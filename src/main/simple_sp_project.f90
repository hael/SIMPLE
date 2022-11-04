module simple_sp_project
include 'simple_lib.f08'
implicit none

public :: sp_project, oritype2segment
private
#include "simple_local_flags.inc"

integer(kind(ENUM_ORISEG)), parameter :: MAXN_OS_SEG = 13

type sp_project
    ! ORIS REPRESENTATIONS OF BINARY FILE SEGMENTS
    ! segments 1-10 reserved for simple program outputs, orientations and files
    ! In segment 7 we stash class averages, ranked class averages, final volumes etc.
    type(oris) :: os_mic    ! micrographs,              segment 1
    type(oris) :: os_stk    ! per-micrograph stack os,  segment 2
    type(oris) :: os_ptcl2D ! per-particle 2D os,       segment 3
    type(oris) :: os_cls2D  ! per-cluster 2D os,        segment 4
    type(oris) :: os_cls3D  ! per-cluster 3D os,        segment 5
    type(oris) :: os_ptcl3D ! per-particle 3D os,       segment 6
    type(oris) :: os_out    ! critical project outputs, segment 7
    type(oris) :: os_optics ! optics groups,            segment 8

    ! ORIS REPRESENTATIONS OF PROJECT DATA / DISTRIBUTED SYSTEM INFO / SYSTEM MANAGEMENT STUFF
    ! segments 11-20 reserved for project info, job management etc.
    type(oris) :: projinfo  ! project information      segment 11
    type(oris) :: jobproc   ! jobid + PID + etc.       segment 12
    type(oris) :: compenv   ! computing environment    segment 13

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
    procedure          :: map_ptcl_ind2stk_ind
    procedure          :: map_cavgs_selection
    ! os_mic related methods
    procedure          :: add_single_movie
    procedure          :: add_movies
    procedure          :: add_intgs
    procedure          :: get_movies_table
    procedure          :: get_mics_table
    procedure          :: get_micparams
    ! os_stk related methods
    procedure          :: add_stk
    procedure, private :: add_stktab_1
    procedure, private :: add_stktab_2
    generic            :: add_stktab => add_stktab_1, add_stktab_2
    procedure          :: add_single_stk
    procedure          :: get_stkname
    procedure          :: get_stkname_and_ind
    procedure, private :: add_scale_tag
    ! os_out related methods
    procedure          :: add_cavgs2os_out
    procedure          :: add_fsc2os_out
    procedure          :: add_frcs2os_out
    procedure          :: add_vol2os_out
    procedure          :: add_entry2os_out
    procedure          :: get_cavgs_stk
    procedure          :: get_vol
    procedure          :: get_fsc
    procedure          :: get_frcs
    procedure          :: get_imginfo_from_osout
    ! getters
    procedure          :: get_n_insegment
    procedure          :: get_n_insegment_state
    procedure          :: get_nptcls
    procedure          :: get_box
    procedure          :: has_boxcoords
    procedure          :: get_boxcoords
    procedure          :: get_smpd
    procedure          :: get_nmovies
    procedure          :: get_nintgs
    procedure          :: get_nframes
    procedure          :: get_nstks
    procedure          :: get_ctfflag
    procedure          :: get_ctfflag_type
    procedure          :: has_phaseplate
    procedure          :: has_boxfile
    procedure          :: get_ctfparams
    procedure          :: get_sp_oris
    procedure          :: ptr2oritype
    procedure          :: is_virgin_field
    procedure          :: get_mic2stk_inds
    ! modifiers
    procedure          :: split_stk
    procedure          :: set_sp_oris
    procedure          :: scale_projfile
    procedure          :: merge_algndocs
    procedure          :: map2Dshifts23D
    procedure          :: map2ptcls
    procedure          :: map2ptcls_state
    procedure          :: replace_project
    procedure          :: merge_stream_projects
    procedure          :: report_state2stk
    procedure          :: set_boxcoords
    ! I/O
    ! printers
    procedure          :: print_info
    procedure          :: print_segment
    ! readers
    procedure          :: read
    procedure          :: read_mic_stk_ptcl2D_segments
    procedure          :: read_non_data_segments
    procedure          :: read_ctfparams_state_eo
    procedure          :: read_segment
    procedure, private :: segreader
    procedure          :: read_segments_info
    procedure          :: read_data_info
    ! writers
    procedure          :: write
    procedure          :: write_segment_inside
    procedure          :: write_non_data_segments
    procedure          :: write_segment2txt
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
                call self%os_mic%new(n,    is_ptcl=.false.)
                os_ptr => self%os_mic
            case('stk')
                call self%os_stk%new(n,    is_ptcl=.false.)
                os_ptr => self%os_stk
            case('ptcl2D')
                call self%os_ptcl2D%new(n, is_ptcl=.true.)
                os_ptr => self%os_ptcl2D
            case('cls2D')
                call self%os_cls2D%new(n,  is_ptcl=.false.)
                os_ptr => self%os_cls2D
            case('cls3D')
                call self%os_cls3D%new(n,  is_ptcl=.false.)
                os_ptr => self%os_cls3D
            case('ptcl3D')
                call self%os_ptcl3D%new(n, is_ptcl=.true.)
                os_ptr => self%os_ptcl3D
            case('optics')
                call self%os_optics%new(n, is_ptcl=.false.)
                os_ptr => self%os_optics
            case('out')
                call self%os_out%new(n,    is_ptcl=.false.)
                os_ptr => self%os_out
            case DEFAULT
                THROW_HARD('unsupported oritype: '//trim(oritype)//'; new_seg_with_ptr')
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
            call self%projinfo%new(1, is_ptcl=.false.)
        endif
        ! projname & profile
        if( self%projinfo%isthere('projname').and.cline%defined('projname') )then
            projname = cline%get_carg('projname')
            call self%projinfo%set(1, 'projname', trim(projname))
            call self%projinfo%set(1, 'projfile', trim(projname)//'.simple')
        else
            if( .not. cline%defined('projname') .and. .not. cline%defined('projfile') )then
                THROW_HARD('the project needs a name, inputted via projname or projfile; update_projinfo')
            endif
            if( cline%defined('projfile') )then
                projfile = cline%get_carg('projfile')
                select case(fname2format(projfile))
                    case('O')
                        call self%projinfo%set(1, 'projfile', trim(projfile) )
                    case DEFAULT
                        THROW_HARD('unsupported format of projfile: '//trim(projfile)//'; update_projinfo')
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
        character(len=:), allocatable    :: projname, qsnam
        integer :: iostat
        if( self%compenv%get_noris() == 1 )then
            ! no need to construct field
        else
            call self%compenv%new(1, is_ptcl=.false.)
        endif
        ! compenv has to be filled as strings as it is used as a string only dictionary
        ! get from environment
        iostat  = simple_getenv('SIMPLE_PATH', env_var)
        if( iostat /= 0 )then
            write(logfhandle,*) 'ERROR! SIMPLE_PATH is not defined in your shell environment!'
            write(logfhandle,*) 'Please refer to installation documentation for correct system configuration'
            stop
        else
            call self%compenv%set(1, 'simple_path', trim(env_var))
        endif
        if( cline%defined('qsys_name') )then
            qsnam = cline%get_carg('qsys_name')
            call self%compenv%set(1, 'qsys_name', trim(qsnam))
            iostat = 0
        else
            iostat = simple_getenv('SIMPLE_QSYS', env_var)
            if( iostat == 0 ) call self%compenv%set(1, 'qsys_name', trim(env_var))
        endif
        if( iostat /= 0 ) THROW_HARD('SIMPLE_QSYS is not defined in your environment; update_compenv')
        iostat = simple_getenv('SIMPLE_EMAIL', env_var)
        if( iostat/=0 ) env_var = 'my.name@uni.edu'
        ! get from command line
        call self%compenv%set(1, 'user_email', trim(env_var))
        if( cline%defined('user_email') )then
            call self%compenv%set(1, 'user_email', cline%get_carg('user_email'))
        else
            call self%compenv%set(1, 'user_email', trim(env_var))
        endif
        if( cline%defined('time_per_image') )then
            call self%compenv%set(1, 'time_per_image', cline%get_rarg('time_per_image'))
        else
            if( .not. self%compenv%isthere('time_per_image') )then
                call self%compenv%set(1, 'time_per_image', real(TIME_PER_IMAGE_DEFAULT))
            endif
        endif
        if( cline%defined('walltime') )then
            call self%compenv%set(1, 'walltime', cline%get_rarg('walltime'))
        else
            if( .not. self%compenv%isthere('walltime') )then
                call self%compenv%set(1, 'walltime', real(WALLTIME_DEFAULT))
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
            call self%compenv%set(1, 'job_memory_per_task', cline%get_rarg('job_memory_per_task') )
        else
            if( .not. self%compenv%isthere('job_memory_per_task') )then
                call self%compenv%set(1, 'job_memory_per_task', real(JOB_MEMORY_PER_TASK_DEFAULT))
            endif
        endif
    end subroutine update_compenv

    !> append segment to current project. BOTH projects must be read in first!
    subroutine append_project( self, proj, oritype )
        class(sp_project), target, intent(inout) :: self, proj
        character(len=*),          intent(in)    :: oritype
        class(oris),          pointer :: os_ptr, os_append_ptr
        type(ori)                     :: o
        type(ctfparams)               :: ctfvar
        character(len=:), allocatable :: stk
        real                          :: smpd, smpd_self, dfx,dfy,angast,phshift
        integer                       :: boxcoords(2),i,iptcl,cnt,n,n2append,nptcls,istate
        select case(trim(oritype))
            case('mic')
                os_ptr => self%os_mic
                os_append_ptr => proj%os_mic
            case('stk')
                os_ptr => self%os_stk
                os_append_ptr => proj%os_stk
            case DEFAULT
                THROW_HARD('oritype: '//trim(oritype)//' unsupported for this purpose; append_project')
        end select
        n2append = os_append_ptr%get_noris()
        if( n2append == 0 )return
        smpd   = os_append_ptr%get(1, 'smpd')
        n      = os_ptr%get_noris()
        if( n == 0 )then
            ! first entry
        else
            smpd_self = os_ptr%get(1, 'smpd')
            if( abs(smpd-smpd_self) > 0.001 )then
                write(logfhandle,*) 'smpd self', smpd_self
                write(logfhandle,*) 'smpd 2 append', smpd
                THROW_HARD('only a project with the same smpd can be appended to the project; append_project')
            endif
        endif
        ! dealing with states
        istate = 1
        if( os_append_ptr%isthere(1,'state') ) istate = os_append_ptr%get_state(1)
        call os_append_ptr%set(1,'state',real(istate))
        ! appending
        select case(trim(oritype))
            case('mic')
                if( n == 0 )then
                    call os_ptr%copy(os_append_ptr, is_ptcl=.false.)
                else
                    ! append
                    call os_ptr%reallocate(n+n2append)
                    cnt = n
                    do i=1,n2append
                        cnt = cnt + 1
                        call os_append_ptr%get_ori(i, o)
                        call os_ptr%set_ori(cnt, o)
                    enddo
                endif
            case('stk')
                ! this assumes there's only one stack in the project to append
                call os_append_ptr%getter(1, 'stk', stk)
                ctfvar = proj%get_ctfparams('ptcl2D', 1)
                cnt    = self%os_ptcl2D%get_noris()
                call self%add_stk(stk, ctfvar)
                nptcls = proj%os_ptcl2D%get_noris()
                do iptcl=1,nptcls
                    cnt = cnt + 1
                    if( proj%has_boxcoords(iptcl) )then
                        call proj%get_boxcoords(iptcl, boxcoords)
                        call self%set_boxcoords(cnt, boxcoords)
                    endif
                    if( ctfvar%ctfflag /= CTFFLAG_NO )then
                        dfx    = proj%os_ptcl2D%get_dfx(iptcl)
                        dfy    = proj%os_ptcl2D%get_dfy(iptcl)
                        angast = proj%os_ptcl2D%get(iptcl, 'angast')
                        call self%os_ptcl2D%set_dfx(cnt,dfx)
                        call self%os_ptcl2D%set_dfy(cnt,dfy)
                        call self%os_ptcl3D%set_dfx(cnt,dfx)
                        call self%os_ptcl3D%set_dfy(cnt,dfy)
                        call self%os_ptcl2D%set(cnt,'angast',angast)
                        call self%os_ptcl3D%set(cnt,'angast',angast)
                        if( proj%os_ptcl2D%isthere(iptcl,'phshift') )then
                            phshift = proj%os_ptcl2D%get(iptcl, 'phshift')
                            call self%os_ptcl2D%set(cnt,'phshift',phshift)
                            call self%os_ptcl3D%set(cnt,'phshift',phshift)
                        endif
                    endif
                enddo
        end select
        nullify(os_ptr, os_append_ptr)
        call o%kill
    end subroutine append_project

    subroutine append_job_descr2jobproc( self, exec_dir, job_descr, did_update )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: exec_dir
        class(chash),      intent(inout) :: job_descr
        logical,           intent(out)   :: did_update
        character(len=:), allocatable    :: edir
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
            call self%jobproc%new(1, is_ptcl=.false.)
            ind = 1
        endif
        call self%jobproc%set_ori(ind,o)
        call o%kill
    end subroutine append_job_descr2jobproc

    ! index management

    subroutine map_ptcl_ind2stk_ind( self, oritype, iptcl, stkind, ind_in_stk )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,                   intent(in)    :: iptcl
        integer,                   intent(out)   :: stkind
        integer,                   intent(out)   :: ind_in_stk
        class(oris), pointer                     :: ptcl_field
        real    :: smpd
        integer :: nptcls, fromp, top, box
        nullify(ptcl_field)
        ! set field pointer
        select case(trim(oritype))
            case('cls3D')
                call self%get_imginfo_from_osout(smpd, box, nptcls)
                if( iptcl < 1 .or. iptcl > nptcls )then
                    write(logfhandle,*) 'iptcl : ', iptcl
                    write(logfhandle,*) 'ncls: ',   nptcls
                    THROW_HARD('iptcl index out of range 1; map_ptcl_ind2stk_ind')
                endif
                stkind     = 1
                ind_in_stk = iptcl
                return
            case('ptcl2D')
                ptcl_field => self%os_ptcl2D
            case('ptcl3D')
                ptcl_field => self%os_ptcl3D
            case DEFAULT
                THROW_HARD('oritype: '//trim(oritype)//' not supported by map_ptcl_ind2stk_ind')
        end select
        nptcls = ptcl_field%get_noris()
        ! first sanity check, range
        if( iptcl < 1 .or. iptcl > nptcls )then
            write(logfhandle,*) 'iptcl : ', iptcl
            write(logfhandle,*) 'nptcls: ', nptcls
            THROW_HARD('iptcl index out of range 2; map_ptcl_ind2stk_ind')
        endif
        ! second sanity check, stack index present in ptcl_field
        if( .not. ptcl_field%isthere(iptcl, 'stkind') )then
            write(logfhandle,*) 'iptcl: ', iptcl
            THROW_HARD('stkind not present in field: '//trim(oritype)//'; map_ptcl_ind2stk_ind')
        endif
        stkind = nint(ptcl_field%get(iptcl, 'stkind'))
        if( ptcl_field%isthere(iptcl, 'indstk') )then
            ind_in_stk = nint(ptcl_field%get(iptcl, 'indstk'))
        else
            ! third sanity check, particle index in range
            fromp = nint(self%os_stk%get(stkind, 'fromp'))
            top   = nint(self%os_stk%get(stkind, 'top'))
            if( iptcl < fromp .or. iptcl > top )then
                write(logfhandle,*) 'iptcl            : ', iptcl
                write(logfhandle,*) 'stkind           : ', stkind
                write(logfhandle,*) 'prange for micstk: ', fromp, top
                THROW_HARD('iptcl index out of micstk range; map_ptcl_ind2stk_ind')
            endif
            ! output index in stack
            ind_in_stk = iptcl - fromp + 1
        endif
        ! cleanup
        nullify(ptcl_field)
    end subroutine map_ptcl_ind2stk_ind

    subroutine map_cavgs_selection( self, states )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: states(:)
        integer, allocatable :: pinds(:)
        integer :: icls, sz_cls2D, sz_states, sz_cls3D, ncls, i
        real    :: rstate
        sz_cls2D  = self%os_cls2D%get_noris()
        sz_states = size(states)
        if( sz_cls2D /= sz_states )then
            write(logfhandle,*) 'size(cls2D): ', sz_cls2D
            write(logfhandle,*) 'sz_states  : ', sz_states
            THROW_HARD('size(cls2D) not consistent with size(states) in map_cavgs_selection, aborting...')
        endif
        ! map selection to self%os_cls2D
        do icls=1,sz_cls2D
            call self%os_cls2D%set(icls, 'state', real(states(icls)))
        end do
        ! map selection to self%os_cls3D
        sz_cls3D = self%os_cls3D%get_noris()
        if( sz_cls3D /= sz_cls2D ) call self%os_cls3D%new(sz_cls2D, is_ptcl=.false.)
        sz_cls3D = sz_cls2D
        do icls=1,sz_cls3D
            call self%os_cls3D%set(icls, 'state', real(states(icls)))
        end do
        ! map selection to self%os_ptcl2D & os_ptcl3D
        ncls = sz_states
        if( self%os_ptcl2D%get_noris() > 0 .and. self%os_ptcl3D%get_noris() > 0)then
            do icls=1,ncls
                call self%os_ptcl2D%get_pinds(icls, 'class', pinds)
                if( allocated(pinds) )then
                    rstate = real(states(icls))
                    do i=1,size(pinds)
                        call self%os_ptcl2D%set(pinds(i), 'state', rstate)
                        call self%os_ptcl3D%set(pinds(i), 'state', rstate)
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
        character(len=LONGSTRLEN)     :: rel_fname
        integer :: n_os_mic, ldim(3), nframes
        ! oris object pointer
        os_ptr => self%os_mic
        ! check that os_mic field is empty
        n_os_mic = os_ptr%get_noris()
        if( n_os_mic > 0 )then
            write(logfhandle,*) 'mic field (self%os_mic) already populated with # entries: ', n_os_mic
            THROW_HARD('add_single_movie')
        endif
        ! update ori
        call os_ptr%new(1, is_ptcl=.false.)
        call make_relativepath(CWD_GLOB, moviename, rel_fname)
        call find_ldim_nptcls(trim(rel_fname), ldim, nframes)
        if( nframes <= 0 )then
            THROW_WARN('# frames in movie: '//trim(rel_fname)//' <= zero, omitting')
        else if( nframes > 1 )then
            call os_ptr%set(1, 'movie', trim(rel_fname))
            call os_ptr%set(1, 'imgkind', 'movie')
            call os_ptr%set(1, 'nframes',    real(nframes))
        else
            call os_ptr%set(1, 'intg',  trim(rel_fname))
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
            case(CTFFLAG_NO)
                call os_ptr%set(1, 'ctf', 'no')
            case(CTFFLAG_YES)
                call os_ptr%set(1, 'ctf', 'yes')
            case(CTFFLAG_FLIP)
                call os_ptr%set(1, 'ctf', 'flip')
            case DEFAULT
                THROW_HARD('ctfflag: '//int2str(ctfvars%ctfflag)//' unsupported; add_single_movie')
        end select
        call os_ptr%set(1,'state',1.) ! default on import
    end subroutine add_single_movie

    !> Add/append movies or micrographs without ctf parameters
    subroutine add_movies( self, movies_array, ctfvars, singleframe )
        class(sp_project), target, intent(inout) :: self
        character(LONGSTRLEN),     intent(in)    :: movies_array(:)
        type(ctfparams),           intent(in)    :: ctfvars
        logical,         optional, intent(in)    :: singleframe
        class(oris),      pointer     :: os_ptr
        type(ctfparams)               :: prev_ctfvars
        character(len=:), allocatable :: name
        character(len=LONGSTRLEN)     :: rel_moviename
        integer                       :: ldim_orig(3), imic, ldim(3), nframes, nmics, nprev_mics, cnt, ntot, nframes_first
        logical                       :: is_movie, l_singleframe
        is_movie = .true.
        l_singleframe = .false.
        if( present(singleframe) ) l_singleframe = singleframe
        ! oris object pointer
        os_ptr => self%os_mic
        ! read movie names
        nmics = size(movies_array)
        ! update oris
        nprev_mics = os_ptr%get_noris()
        ntot       = nmics + nprev_mics
        if( nprev_mics == 0 )then
            call os_ptr%new(ntot, is_ptcl=.false.)
        else
            prev_ctfvars = self%get_micparams(1)
            if( ctfvars%ctfflag /= prev_ctfvars%ctfflag ) THROW_HARD('CTF infos do not match! add_movies')
            if( .not.is_equal(ctfvars%smpd, prev_ctfvars%smpd )) THROW_HARD('The sampling distances do not match! add_movies')
            if( .not.is_equal(ctfvars%cs,   prev_ctfvars%cs   )) THROW_HARD('The spherical aberrations do not match! add_movies')
            if( .not.is_equal(ctfvars%kv,   prev_ctfvars%kv   )) THROW_HARD('The voltages do not match! add_movies')
            if( .not.is_equal(ctfvars%fraca,prev_ctfvars%fraca)) THROW_HARD('The amplitude contrasts do not match! add_movies')
            if( ctfvars%l_phaseplate.neqv.prev_ctfvars%l_phaseplate ) THROW_HARD('Phaseplate infos do not match! add_movies')
            call os_ptr%reallocate(ntot)
        endif
        cnt = 0
        nframes_first = 0
        do imic=nprev_mics + 1,ntot
            cnt = cnt + 1
            call make_relativepath(CWD_GLOB,movies_array(cnt),rel_moviename)
            call find_ldim_nptcls(trim(rel_moviename), ldim, nframes)
            if( cnt == 1 )then
                ldim_orig = ldim
            else
                if( ldim(1) /= ldim_orig(1) .or. ldim(2) /= ldim_orig(2) )then
                    write(logfhandle,*)'Inconsistent size for file: ',trim(movies_array(cnt))
                    write(logfhandle,*)'Dimensions: ', ldim(1),'x ',ldim(2), ' vs. previous dimensions: ', ldim_orig(1),'x ',ldim_orig(2)
                    THROW_HARD('All files imported must have identical diemnsions!')
                endif
            endif
            if( nframes <= 0 )then
                THROW_WARN('# frames in movie: '//trim(movies_array(imic))//' <= zero, omitting')
                cycle
            else
                if( nframes > 1 )then
                    call os_ptr%set(imic, 'movie', trim(rel_moviename))
                    call os_ptr%set(imic, 'imgkind', 'movie')
                    is_movie = .true.
                else
                    if( l_singleframe )then
                        call os_ptr%set(imic, 'frame',   trim(rel_moviename))
                        call os_ptr%set(imic, 'imgkind', 'frame')
                    else
                        call os_ptr%set(imic, 'intg',    trim(rel_moviename))
                        call os_ptr%set(imic, 'imgkind', 'mic')
                    endif
                    is_movie = .false.
                endif
                if( nframes_first == 0 )then
                    nframes_first = nframes
                else
                    if( nframes /= nframes_first )then
                        write(logfhandle,*) trim(rel_moviename), ' has ', nframes, ' frame(s)'
                        write(logfhandle,*) 'Previous import have ', nframes_first, ' frame(s)'
                        THROW_HARD('You cannot import both micrographs and movies at the same time! add_movies')
                    endif
                endif
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
                case(CTFFLAG_NO)
                    call os_ptr%set(imic, 'ctf', 'no')
                case(CTFFLAG_YES)
                    call os_ptr%set(imic, 'ctf', 'yes')
                case(CTFFLAG_FLIP)
                    call os_ptr%set(imic, 'ctf', 'flip')
                case DEFAULT
                    THROW_HARD('unsupported ctfflag: '//int2str(ctfvars%ctfflag)//'; add_movies')
            end select
        enddo
        if( is_movie )then
            name = 'MOVIE(S)'
        else
            if( l_singleframe )then
                name = 'FRAME(S)'
            else
                name = 'MICROGRAPH(S)'
            endif
        endif
        write(logfhandle,'(A13,I6,A1,A)')'>>> IMPORTED ', nmics,' ', trim(name)
        write(logfhandle,'(A20,A,A1,I6)')'>>> TOTAL NUMBER OF ', trim(name),':',ntot
    end subroutine add_movies

    !> Add/append micrographs with ctf parameters
    subroutine add_intgs( self, intgs_array, os, ctfvars )
        class(sp_project), target, intent(inout) :: self
        character(LONGSTRLEN),     intent(in)    :: intgs_array(:)
        class(oris),               intent(in)    :: os
        type(ctfparams),           intent(in)    :: ctfvars
        type(ctfparams)           :: prev_ctfvars, ctfparms
        character(len=LONGSTRLEN) :: rel_micname
        real                      :: intg_smpd
        integer                   :: imic,ldim(3),nframes,nintgs,nprev_intgs,nprev_mics,cnt,ntot
        nprev_mics  = self%os_mic%get_noris()
        nprev_intgs = self%get_nintgs()
        if( nprev_mics > 0 )then
            if( nprev_mics /= nprev_intgs )then
                THROW_HARD('Cannot add lone micrographs to a project with movies; add_intgs')
            endif
            if( nprev_intgs == 0 )then
                THROW_HARD('Cannot add micrographs to a project with movies only; add_intgs')
            endif
            ! previous micrographs parameters
            prev_ctfvars = self%os_mic%get_ctfvars(1)
            if(.not.is_equal(ctfvars%smpd, prev_ctfvars%smpd )) THROW_HARD('Inconsistent sampling distance; add_intgs')
            if(.not.is_equal(ctfvars%cs,   prev_ctfvars%cs   )) THROW_HARD('Inconsistent spherical aberration; add_intgs')
            if(.not.is_equal(ctfvars%kv,   prev_ctfvars%kv   )) THROW_HARD('Inconsistent voltage; add_intgs')
            if(.not.is_equal(ctfvars%fraca,prev_ctfvars%fraca)) THROW_HARD('Inconsistent amplituce contrast; add_intgs')
            if(ctfvars%ctfflag /= prev_ctfvars%ctfflag) THROW_HARD('Incompatible CTF flag; add_intgs')
            if(ctfvars%l_phaseplate .neqv. prev_ctfvars%l_phaseplate ) THROW_HARD('Incompatible phaseplate info; add_intgs')
        endif
        ! read movie names
        nintgs = size(intgs_array)
        if( nintgs /= os%get_noris() )then
            THROW_HARD('Inconsistent # of mics & ctf parameters; add_intgs')
        endif
        ! update oris
        if( nprev_intgs == 0 )then
            ! first import
            call self%os_mic%new(nintgs, is_ptcl=.false.)
            ntot = nintgs
        else
            ! append
            ntot = nintgs+nprev_intgs
            call self%os_mic%reallocate(ntot)
        endif
        cnt = 0
        do imic=nprev_intgs+1,ntot
            cnt = cnt + 1
            call make_relativepath(CWD_GLOB,intgs_array(cnt),rel_micname)
            call find_ldim_nptcls(trim(rel_micname), ldim, nframes, smpd=intg_smpd)
            if( nframes <= 0 )then
                THROW_HARD('# frames in movie: '//trim(intgs_array(cnt))//' <= zero; add_intgs')
            else if( nframes > 1 )then
                THROW_HARD('Not the interface for adding movies; add_intgs')
            endif
            if( nprev_intgs > 0 )then
                if( .not.is_equal(intg_smpd,prev_ctfvars%smpd) )then
                    THROW_HARD('Incompatible sampling distance: '//trim(intgs_array(cnt))//'; add_intgs')
                endif
            endif
            ! updates segment
            ctfparms = os%get_ctfvars(cnt)
            call self%os_mic%set(imic, 'intg', rel_micname)
            call self%os_mic%set(imic, 'imgkind', 'mic')
            call self%os_mic%set(imic, 'xdim',    real(ldim(1)))
            call self%os_mic%set(imic, 'ydim',    real(ldim(2)))
            call self%os_mic%set(imic, 'smpd',    ctfvars%smpd)
            call self%os_mic%set(imic, 'kv',      ctfvars%kv)
            call self%os_mic%set(imic, 'cs',      ctfvars%cs)
            call self%os_mic%set(imic, 'fraca',   ctfvars%fraca)
            if( os%isthere(cnt,'state') )then
                call self%os_mic%set(imic, 'state', os%get(cnt,'state'))
            else
                call self%os_mic%set(imic, 'state', 1.0)
            endif
            if( ctfvars%l_phaseplate )then
                call self%os_mic%set(imic, 'phaseplate', 'yes')
            else
                call self%os_mic%set(imic, 'phaseplate', 'no')
            endif
            select case(ctfvars%ctfflag)
                case(CTFFLAG_NO)
                    call self%os_mic%set(imic, 'ctf', 'no')
                case(CTFFLAG_YES)
                    call self%os_mic%set(imic, 'ctf',    'yes')
                    call self%os_mic%set_dfx(imic,       ctfparms%dfx)
                    call self%os_mic%set_dfy(imic,       ctfparms%dfy)
                    call self%os_mic%set(imic, 'angast', ctfparms%angast)
                    call self%os_mic%set(imic, 'phshift',ctfparms%phshift)
                case(CTFFLAG_FLIP)
                    call self%os_mic%set(imic, 'ctf', 'flip')
            end select
        enddo
        write(logfhandle,'(A,I6,A)')'>>> IMPORTED ', nintgs,' INTEGRATED MOVIES'
        write(logfhandle,'(A,I6)')'>>> TOTAL NUMBER OF MICROGRAPHS:',ntot
    end subroutine add_intgs

    ! returns list of movies regardless of 'imgkind' key as it overrides movies
    subroutine get_movies_table( self, moviestab )
        class(sp_project),                      intent(inout) :: self
        character(len=LONGSTRLEN), allocatable, intent(out)   :: moviestab(:)
        character(len=:), allocatable :: mov
        integer :: i,n,cnt
        if(allocated(moviestab))deallocate(moviestab)
        n = 0
        do i=1,self%os_mic%get_noris()
            if(self%os_mic%isthere(i,'movie')) n = n+1
        enddo
        if( n==0 )return
        allocate(moviestab(n))
        cnt = 0
        do i=1,self%os_mic%get_noris()
            if(self%os_mic%isthere(i,'movie'))then
                cnt = cnt + 1
                call self%os_mic%getter(i,'movie',mov)
                moviestab(cnt) = trim(mov)
            endif
        enddo
    end subroutine get_movies_table

    subroutine get_mics_table( self, micstab )
        class(sp_project),                      intent(inout) :: self
        character(len=LONGSTRLEN), allocatable, intent(out)   :: micstab(:)
        character(len=:), allocatable :: imgkind, mic
        integer :: i,n,cnt
        if(allocated(micstab))deallocate(micstab)
        n = self%get_nintgs()
        if( n==0 )return
        allocate(micstab(n))
        cnt = 0
        do i=1,self%os_mic%get_noris()
            if(self%os_mic%isthere(i,'imgkind'))then
                call self%os_mic%getter(i,'imgkind',imgkind)
                if( trim(imgkind).eq.'mic' )then
                    cnt = cnt + 1
                    call self%os_mic%getter(i,'intg',mic)
                    micstab(cnt) = trim(mic)
                endif
            endif
        enddo
    end subroutine get_mics_table

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
        character(len=LONGSTRLEN)     :: stk_relpath, cwd
        integer :: ldim(3), nptcls, n_os_stk, n_os_ptcl2D, n_os_ptcl3D
        integer :: i, fromp, top
        ! full path and existence check
        call simple_getcwd(cwd)
        call make_relativepath(cwd, stk, stk_relpath)
        ! find dimension of inputted stack
        call find_ldim_nptcls(stk_relpath, ldim, nptcls)
        if( ldim(1) /= ldim(2) )then
            write(logfhandle,*) 'xdim: ', ldim(1)
            write(logfhandle,*) 'ydim: ', ldim(2)
            THROW_HARD('nonsquare particle images not supported; add_stk')
        endif
        ! updates_fields
        n_os_stk    = self%os_stk%get_noris() + 1
        n_os_ptcl2D = self%os_ptcl2D%get_noris()
        n_os_ptcl3D = self%os_ptcl3D%get_noris()
        if( n_os_stk == 1 )then
            call self%os_stk%new(1,         is_ptcl=.false.)
            call self%os_ptcl2D%new(nptcls, is_ptcl=.true.)
            call self%os_ptcl3D%new(nptcls, is_ptcl=.true.)
            fromp = 1
            top   = nptcls
        else
            ! stk
            if( .not.self%os_stk%isthere(n_os_stk-1,'top') )then
                THROW_HARD('FROMP/TOP keys should always be informed; add_stk')
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
        call self%os_stk%set(n_os_stk, 'stk',     trim(stk_relpath))
        call self%os_stk%set(n_os_stk, 'box',     real(ldim(1)))
        call self%os_stk%set(n_os_stk, 'nptcls',  real(nptcls))
        call self%os_stk%set(n_os_stk, 'fromp',   real(fromp))
        call self%os_stk%set(n_os_stk, 'top',     real(top))
        call self%os_stk%set(n_os_stk, 'stkkind', 'split')
        call self%os_stk%set(n_os_stk, 'imgkind', 'ptcl')
        call self%os_stk%set(n_os_stk, 'state',   1.0) ! default on import
        select case(ctfvars%ctfflag)
            case(CTFFLAG_NO,CTFFLAG_YES,CTFFLAG_FLIP)
                call self%os_stk%set_ctfvars(n_os_stk, ctfvars)
            case DEFAULT
                THROW_HARD('unsupported ctfflag: '//int2str(ctfvars%ctfflag)//'; add_stk')
        end select
        call self%os_stk%set_ctfvars(n_os_stk, ctfvars)
        ! update particle oris objects
        do i = 1, nptcls
            call o%new(is_ptcl=.true.)
            call o%set_dfx(      ctfvars%dfx)
            call o%set_dfy(      ctfvars%dfy)
            call o%set('angast', ctfvars%angast)
            if( ctfvars%l_phaseplate ) call o%set('phshift', ctfvars%phshift)
            call o%set('stkind', real(n_os_stk))
            call o%set('state',1.) ! default on import
            call self%os_ptcl2D%set_ori(n_os_ptcl2D+i, o)
            call self%os_ptcl3D%set_ori(n_os_ptcl3D+i, o)
        enddo
        call o%kill
    end subroutine add_stk

    subroutine add_single_stk( self, stk, ctfvars, os )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: stk
        type(ctfparams),   intent(in)    :: ctfvars ! CTF parameters associated with stk (smpd,kv,cs,fraca,phaseplate)
        class(oris),       intent(inout) :: os      ! parameters associated with stk (dfx,dfy,angast,phshift)
        character(len=:), allocatable :: stk_abspath, projname, fbody
        character(len=LONGSTRLEN)     :: stk_relpath
        integer                       :: n_os_stk, n_os_ptcl2D, n_os_ptcl3D, ldim(3), nptcls
        call self%projinfo%getter(1, 'projname', projname)
        if( str_has_substr(stk, 'mrc') )then
            fbody = get_fbody(basename(stk), 'mrc')
        else if( str_has_substr(stk, 'mrcs') )then
            fbody = get_fbody(basename(stk), 'mrcs')
        else
            THROW_HARD('Unsupported stack format; use *.mrc or *.mrcs for import')
        endif
        if( str_has_substr(trim(projname), fbody) ) THROW_HARD('stack for import('//trim(stk)//') not allowed to have same name as project')
        ! check that stk field is empty
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk > 0 )then
            write(logfhandle,*) 'stack field (self%os_stk) already populated with # entries: ', n_os_stk
            THROW_HARD('add_single_stk')
        endif
        ! check that particle fields are empty
        n_os_ptcl2D = self%os_ptcl2D%get_noris()
        if( n_os_ptcl2D > 0 )then
            write(logfhandle,*) 'ptcl2D field (self%os_ptcl2D) already populated with # entries: ', n_os_ptcl2D
            THROW_HARD('empty particle fields in project file assumed; add_single_stk')
        endif
        n_os_ptcl3D = self%os_ptcl3D%get_noris()
        if( n_os_ptcl3D > 0 )then
            write(logfhandle,*) 'ptcl3D field (self%os_ptcl3D) already populated with # entries: ', n_os_ptcl3D
            THROW_HARD('empty particle fields in project file assumed; add_single_stk')
        endif
        call self%os_ptcl2D%copy(os, is_ptcl=.true.)
        call self%os_ptcl3D%copy(os, is_ptcl=.true.)
        call self%os_ptcl2D%set_all2single('stkind', 1.)
        call self%os_ptcl3D%set_all2single('stkind', 1.)
        if( .not. self%os_ptcl2D%isthere('state') ) call self%os_ptcl2D%set_all2single('state',  1.) ! default on import
        if( .not. self%os_ptcl3D%isthere('state') ) call self%os_ptcl3D%set_all2single('state',  1.) ! default on import
        ! full path and existence check
        stk_abspath = simple_abspath(stk,'sp_project :: add_single_stk')
        ! find dimension of inputted stack
        call find_ldim_nptcls(trim(stk_abspath), ldim, nptcls)
        if( ldim(1) /= ldim(2) )then
            write(logfhandle,*) 'xdim: ', ldim(1)
            write(logfhandle,*) 'ydim: ', ldim(2)
            THROW_HARD('nonsquare particle images not supported; add_single_stk')
        endif
        ! path
        call make_relativepath(CWD_GLOB, stk_abspath, stk_relpath)
        ! records
        call self%os_stk%new(1, is_ptcl=.false.)
        call self%os_stk%set(1, 'stk',     trim(stk_relpath))
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
            if( .not. os%isthere('phshift') ) THROW_HARD('phaseplate=yes & input oris lack phshift; add_single_stk')
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
                THROW_HARD('unsupported ctfflag: '//int2str(ctfvars%ctfflag)//'; add_single_stk')
        end select
    end subroutine add_single_stk

    ! adds stktab given per-stk parameters
    subroutine add_stktab_1( self, stkfnames, os )
        class(sp_project),     intent(inout) :: self
        character(LONGSTRLEN), intent(inout) :: stkfnames(:)
        class(oris),           intent(inout) :: os ! parameters associated with stktab
        integer,  allocatable :: nptcls_arr(:)
        type(ori)             :: o_ptcl, o_stk
        character(LONGSTRLEN) :: stk_relpath
        integer   :: ldim(3), ldim_here(3), n_os_ptcl2D, n_os_ptcl3D, n_os_stk, istate
        integer   :: i, istk, fromp, top, nptcls, n_os, nstks, nptcls_tot, stk_ind
        nstks = size(stkfnames)
        ! check that inputs are of conforming sizes
        n_os = os%get_noris()
        if( n_os /= nstks )then
            write(logfhandle,*) '# input oris      : ', n_os
            write(logfhandle,*) '# stacks in stktab: ', nstks
            THROW_HARD('nonconforming sizes of inputs; add_stktab')
        endif
        ! first pass for sanity check and determining dimensions
        allocate(nptcls_arr(nstks),source=0)
        do istk=1,nstks
            ! full path and existence check
            call make_relativepath(CWD_GLOB,stkfnames(istk),stk_relpath)
            stkfnames(istk) = trim(stk_relpath)
            call os%get_ori(istk, o_stk)
            ! logical dimension management
            call find_ldim_nptcls(stkfnames(istk), ldim, nptcls)
            ldim(3) = 1
            if( istk == 1 )then
                ldim_here = ldim
            else
                if( .not. all(ldim_here == ldim) )then
                    write(logfhandle,*) 'micrograph stack #  : ', istk
                    write(logfhandle,*) 'stk name            : ', trim(stkfnames(istk))
                    write(logfhandle,*) 'ldim in object      : ', ldim_here
                    write(logfhandle,*) 'ldim read from stack: ', ldim
                    THROW_HARD('inconsistent logical dimensions; add_stktab')
                endif
            endif
            if( ldim(1) /= ldim(2) )then
                write(logfhandle,*) 'stk name: ', trim(stkfnames(istk))
                write(logfhandle,*) 'xdim:     ', ldim(1)
                write(logfhandle,*) 'ydim:     ', ldim(2)
                THROW_HARD('nonsquare particle images not supported; add_stktab')
            endif
            ! check variable presence
            if( .not. o_stk%isthere('ctf') )    THROW_HARD('ERROR! ctf flag missing in os input; add_stktab')
            if( .not. o_stk%isthere('smpd') )   THROW_HARD('ERROR! smpd missing in os input; add_stktab')
            if( .not. o_stk%isthere('kv') )     THROW_HARD('ERROR! kv missing in os input; add_stktab')
            if( .not. o_stk%isthere('cs') )     THROW_HARD('ERROR! cs missing in os input; add_stktab')
            if( .not. o_stk%isthere('fraca') )  THROW_HARD('ERROR! fraca missing in os input; add_stktab')
            if( .not. o_stk%isthere('dfx') )    THROW_HARD('ERROR! dfx missing in os input; add_stktab')
            if( .not. o_stk%isthere('dfy') )    THROW_HARD('ERROR! dfy missing in os input; add_stktab')
            if( .not. o_stk%isthere('angast') ) THROW_HARD('ERROR! angast missing in os input; add_stktab')
            ! stash number of images
            nptcls_arr(istk) = nptcls
        enddo
        ! oris allocation
        nptcls_tot  = sum(nptcls_arr)
        n_os_stk    = self%os_stk%get_noris()
        n_os_ptcl2D = self%os_ptcl2D%get_noris()
        n_os_ptcl3D = self%os_ptcl3D%get_noris()
        if( n_os_stk == 0 )then
            call self%os_stk%new(nstks,         is_ptcl=.false.)
            call self%os_ptcl2D%new(nptcls_tot, is_ptcl=.true.)
            call self%os_ptcl3D%new(nptcls_tot, is_ptcl=.true.)
            fromp = 1
        else
            if( .not.self%os_stk%isthere(n_os_stk,'top') )then
                THROW_HARD('FROMP/TOP keys should always be informed; add_stk')
            endif
            call self%os_stk%reallocate(n_os_stk + nstks)
            call self%os_ptcl2D%reallocate(n_os_ptcl2D + nptcls_tot)
            call self%os_ptcl3D%reallocate(n_os_ptcl3D + nptcls_tot)
            fromp = nint(self%os_stk%get(n_os_stk,'top')) + 1
        endif
        ! parameters transfer
        do istk=1,nstks
            call os%get_ori(istk, o_stk)
            top     = fromp + nptcls_arr(istk) - 1 ! global index
            stk_ind = n_os_stk + istk
            ! updates stk segment
            call self%os_stk%set_ori(stk_ind, o_stk)
            call self%os_stk%set(stk_ind, 'stk',     trim(stkfnames(istk)))
            call self%os_stk%set(stk_ind, 'box',     real(ldim(1)))
            call self%os_stk%set(stk_ind, 'nptcls',  real(nptcls_arr(istk)))
            call self%os_stk%set(stk_ind, 'fromp',   real(fromp))
            call self%os_stk%set(stk_ind, 'top',     real(top))
            call self%os_stk%set(stk_ind, 'stkkind', 'split')
            call self%os_stk%set(stk_ind, 'imgkind', 'ptcl')
            istate = 1 ! default on import
            if( o_stk%isthere('state') ) istate = o_stk%get_state()
            call self%os_stk%set(stk_ind, 'state', real(istate))
            ! updates particles segment
            call o_ptcl%new(is_ptcl=.true.)
            call o_ptcl%set_dfx(      o_stk%get_dfx())
            call o_ptcl%set_dfy(      o_stk%get_dfy())
            call o_ptcl%set('angast', o_stk%get('angast'))
            if( o_stk%isthere('phshift') ) call o_ptcl%set('phshift', o_stk%get('phshift'))
            call o_ptcl%set('stkind', real(stk_ind))
            call o_ptcl%set('state',  real(istate))
            do i=1,nptcls_arr(istk)
                call self%os_ptcl2D%set_ori(fromp+i-1, o_ptcl)
                call self%os_ptcl3D%set_ori(fromp+i-1, o_ptcl)
            enddo
            ! update
            fromp = top + 1 ! global index
        enddo
        call o_ptcl%kill
        call o_stk%kill
    end subroutine add_stktab_1

    ! adds stktab given per-particle parameters
    subroutine add_stktab_2( self, stkfnames, ctfvars, os )
        class(sp_project),     intent(inout) :: self
        character(LONGSTRLEN), intent(inout) :: stkfnames(:)
        type(ctfparams),       intent(in)    :: ctfvars
        class(oris),           intent(inout) :: os ! parameters associated with stktab
        integer,  allocatable :: nptcls_arr(:)
        type(ori)             :: o_ptcl, o_stk
        character(LONGSTRLEN) :: stk_relpath
        integer   :: ldim(3), ldim_here(3), n_os_ptcl2D, n_os_ptcl3D, n_os_stk, istate
        integer   :: i, istk, fromp, top, nptcls, n_os, nstks, nptcls_tot, stk_ind, pind
        nstks = size(stkfnames)
        n_os  = os%get_noris()
        ! first pass for sanity check and determining dimensions
        allocate(nptcls_arr(nstks),source=0)
        do istk=1,nstks
            ! full path and existence check
            call make_relativepath(CWD_GLOB,stkfnames(istk),stk_relpath)
            stkfnames(istk) = trim(stk_relpath)
            ! logical dimension management
            call find_ldim_nptcls(stkfnames(istk), ldim, nptcls)
            ldim(3) = 1
            if( istk == 1 )then
                ldim_here = ldim
            else
                if( .not. all(ldim_here == ldim) )then
                    write(logfhandle,*) 'micrograph stack #  : ', istk
                    write(logfhandle,*) 'stk name            : ', trim(stkfnames(istk))
                    write(logfhandle,*) 'ldim in object      : ', ldim_here
                    write(logfhandle,*) 'ldim read from stack: ', ldim
                    THROW_HARD('inconsistent logical dimensions; add_stktab')
                endif
            endif
            if( ldim(1) /= ldim(2) )then
                write(logfhandle,*) 'stk name: ', trim(stkfnames(istk))
                write(logfhandle,*) 'xdim:     ', ldim(1)
                write(logfhandle,*) 'ydim:     ', ldim(2)
                THROW_HARD('nonsquare particle images not supported; add_stktab')
            endif
            ! stash number of images
            nptcls_arr(istk) = nptcls
        enddo
        nptcls_tot  = sum(nptcls_arr)
        if( n_os /= nptcls_tot )then
            write(logfhandle,*) '# input oris               : ', n_os
            write(logfhandle,*) '# ptcls in stacks in stktab: ', nptcls_tot
            THROW_HARD('nonconforming sizes of inputs; add_stktab_2')
        endif
        ! oris allocation
        n_os_stk    = self%os_stk%get_noris()
        n_os_ptcl2D = self%os_ptcl2D%get_noris()
        n_os_ptcl3D = self%os_ptcl3D%get_noris()
        if( n_os_stk == 0 )then
            call self%os_stk%new(nstks,         is_ptcl=.false.)
            call self%os_ptcl2D%new(nptcls_tot, is_ptcl=.true.)
            call self%os_ptcl3D%new(nptcls_tot, is_ptcl=.true.)
            fromp = 1
        else
            if( .not.self%os_stk%isthere(n_os_stk,'top') )then
                THROW_HARD('FROMP/TOP keys should always be informed; add_stktab_2')
            endif
            call self%os_stk%reallocate(n_os_stk + nstks)
            call self%os_ptcl2D%reallocate(n_os_ptcl2D + nptcls_tot)
            call self%os_ptcl3D%reallocate(n_os_ptcl3D + nptcls_tot)
            fromp = nint(self%os_stk%get(n_os_stk,'top')) + 1
        endif
        ! parameters transfer
        do istk=1,nstks
            top     = fromp + nptcls_arr(istk) - 1 ! global index
            stk_ind = n_os_stk + istk
            ! updates stk segment
            call o_stk%new(is_ptcl=.false.)
            call o_stk%set('stk',     trim(stkfnames(istk)))
            call o_stk%set('box',     real(ldim(1)))
            call o_stk%set('nptcls',  real(nptcls_arr(istk)))
            call o_stk%set('fromp',   real(fromp))
            call o_stk%set('top',     real(top))
            call o_stk%set('stkkind', 'split')
            call o_stk%set('imgkind', 'ptcl')
            call o_stk%set('smpd',    ctfvars%smpd)
            call o_stk%set('kv',      ctfvars%kv)
            call o_stk%set('cs',      ctfvars%cs)
            call o_stk%set('fraca',   ctfvars%fraca)
            call o_stk%set('state',   1.0) ! default on import
            if( ctfvars%l_phaseplate )then
                call o_stk%set('phaseplate', 'yes')
                call o_stk%set('phshift',    ctfvars%phshift)
            else
                call o_stk%set('phaseplate', 'no')
            endif
            select case(ctfvars%ctfflag)
                case(CTFFLAG_NO)
                    call o_stk%set('ctf', 'no')
                case(CTFFLAG_YES)
                    call o_stk%set('ctf', 'yes')
                case(CTFFLAG_FLIP)
                    call o_stk%set('ctf', 'flip')
                case DEFAULT
                    THROW_HARD('unsupported ctfflag: '//int2str(ctfvars%ctfflag)//'; add_stktab_2')
            end select
            call self%os_stk%set_ori(stk_ind, o_stk)
            ! updates particles segment
            do i=1,nptcls_arr(istk)
                pind = fromp+i-1
                call os%get_ori(pind, o_ptcl)
                call o_ptcl%set('stkind', real(stk_ind))
                istate = 1
                if( o_ptcl%isthere('state') ) istate = o_ptcl%get_state()
                call o_ptcl%set('state',  real(istate))
                select case(ctfvars%ctfflag)
                    case(CTFFLAG_YES,CTFFLAG_FLIP)
                        if( .not.o_ptcl%isthere('dfx') )then
                            call o_ptcl%print_ori
                            THROW_HARD('Missing defocus parameter(s) for particle: '//int2str(pind))
                        endif
                    case DEFAULT
                        ! all good
                end select
                if( ctfvars%l_phaseplate )then
                    if( .not.o_ptcl%isthere('phshift') )then
                        call o_ptcl%print_ori
                        THROW_HARD('Missing phase-shift parameter for particle: '//int2str(pind))
                    endif
                endif
                call self%os_ptcl2D%set_ori(pind, o_ptcl)
                call self%os_ptcl3D%set_ori(pind, o_ptcl)
            enddo
            ! update
            fromp = top + 1 ! global index
        enddo
        call o_ptcl%kill
        call o_stk%kill
    end subroutine add_stktab_2

    !>  Only commits to disk when a change to the project is made
    subroutine split_stk( self, nparts, dir )
        use simple_map_reduce, only: split_nobjs_even
        use simple_image,      only: image
        use simple_stack_io,   only: stack_io
        class(sp_project),          intent(inout) :: self
        integer,                    intent(in)    :: nparts
        character(len=*), optional, intent(in)    :: dir
        character(len=*), parameter   :: EXT = '.mrc'
        type(image)                   :: img
        type(ori)                     :: orig_stk
        type(stack_io)                :: stkio_r, stkio_w
        character(len=:), allocatable :: stk, tmp_dir, stkpart
        character(len=:), allocatable :: dest_stkpart
        character(len=LONGSTRLEN) :: stk_relpath, cwd
        integer :: parts(nparts,2), ind_in_stk, iptcl, cnt, istk, box, n_os_stk
        integer :: nptcls, nptcls_part, numlen
        real    :: smpd
        if( nparts < 2 )return
        ! check that stk field is not empty
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk == 0 )then
            THROW_HARD('No stack to split! split_stk')
        else if( n_os_stk > 1 )then ! re-splitting not supported
            return
        endif
        ! get original simple_parameters
        call self%os_stk%get_ori(1, orig_stk)
        ! copy prep
        nptcls  = self%get_nptcls()
        parts   = split_nobjs_even( nptcls, nparts )
        numlen  = len_trim(int2str(nparts))
        ! images copy
        smpd = orig_stk%get('smpd')
        box  = nint(orig_stk%get('box'))
        call img%new([box,box,1], smpd)
        call simple_getcwd(cwd)
        if( present(dir) )then
            tmp_dir = filepath(trim(dir),'tmp_stacks')
        else
            tmp_dir = filepath(trim(cwd),'tmp_stacks')
        endif
        call simple_mkdir(trim(tmp_dir),errmsg="sp_project::split_stk")
        write(logfhandle,'(a)') '>>> SPLITTING STACK INTO PARTS'
        ! just to get the name of the stack to read from
        call self%get_stkname_and_ind('ptcl2D', 1, stk, ind_in_stk)
        call stkio_r%open(stk, smpd, 'read')
        do istk = 1,nparts
            call progress(istk,nparts)
            stkpart = filepath(trim(tmp_dir),'stack_part'//int2str_pad(istk,numlen)//EXT)
            call stkio_w%open(stkpart, smpd, 'write', box=box, is_ft=.false.)
            cnt = 0
            do iptcl = parts(istk,1), parts(istk,2)
                cnt = cnt + 1
                call self%get_stkname_and_ind('ptcl2D', iptcl, stk, ind_in_stk)
                call stkio_r%read(ind_in_stk, img)
                call stkio_w%write(cnt, img)
            enddo
            deallocate(stkpart)
            call stkio_w%close
        enddo
        call stkio_r%close
        call img%kill
        call self%os_stk%new(nparts, is_ptcl=.false.)
        if( present(dir) )then
           call simple_mkdir(filepath(trim(dir),trim(STKPARTSDIR)),errmsg="sp_project::split_stk")
        else
           call simple_mkdir(trim(STKPARTSDIR),errmsg="sp_project::split_stk")
        endif
        do istk = 1,nparts
            ! file stuff
            stkpart = filepath(trim(tmp_dir),'stack_part'//int2str_pad(istk,numlen)//EXT)
            if( present(dir) )then
                dest_stkpart = filepath(trim(dir),trim(STKPARTFBODY)//int2str_pad(istk,numlen)//EXT)
            else
                allocate(dest_stkpart, source=trim(STKPARTFBODY)//int2str_pad(istk,numlen)//EXT)
            endif
            call simple_rename(trim(stkpart), trim(dest_stkpart))
            call make_relativepath(cwd, dest_stkpart, stk_relpath)
            nptcls_part = parts(istk,2) - parts(istk,1) + 1
            ! set original before overriding
            call self%os_stk%set_ori(istk, orig_stk)
            ! override
            call self%os_stk%set(istk, 'stk',     trim(stk_relpath))
            call self%os_stk%set(istk, 'nptcls',  real(nptcls_part))
            call self%os_stk%set(istk, 'fromp',   real(parts(istk,1)))
            call self%os_stk%set(istk, 'top',     real(parts(istk,2)))
            call self%os_stk%set(istk, 'stkkind', 'split')
            do iptcl=parts(istk,1),parts(istk,2)
                call self%os_ptcl2D%set(iptcl,'stkind',real(istk))
                call self%os_ptcl3D%set(iptcl,'stkind',real(istk))
            enddo
            deallocate(stkpart, dest_stkpart)
        enddo
        call self%write
        call simple_rmdir(tmp_dir,errmsg="sp_project::split_stk")
        call orig_stk%kill
    end subroutine split_stk

    function get_stkname( self, imic ) result( stkname )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: imic
        character(len=:), allocatable    :: stkname
        integer :: nmics
        nmics = self%os_stk%get_noris()
        if( imic < 1 .or. imic > nmics )then
            write(logfhandle,*) 'imic : ', imic
            write(logfhandle,*) 'nmics: ', nmics
            THROW_HARD('imic index out of range; get_stkname')
        endif
        stkname = trim(self%os_stk%get_static(imic, 'stk'))
    end function get_stkname

    subroutine get_stkname_and_ind( self, oritype, iptcl, stkname, ind_in_stk )
        class(sp_project), target,     intent(inout) :: self
        character(len=*),              intent(in)    :: oritype
        integer,                       intent(in)    :: iptcl
        character(len=:), allocatable, intent(out)   :: stkname
        integer,                       intent(out)   :: ind_in_stk
        real    :: smpd
        integer :: stkind, ncls
        ! do the index mapping
        call self%map_ptcl_ind2stk_ind(oritype, iptcl, stkind, ind_in_stk )
        ! output name
        if( allocated(stkname) ) deallocate(stkname)
        if( trim(oritype) .eq. 'cls3D' ) then
            call self%get_cavgs_stk(stkname, ncls, smpd)
        else
            stkname = trim(self%os_stk%get_static(stkind, 'stk'))
        endif
    end subroutine get_stkname_and_ind

    subroutine add_scale_tag( self, dir )
        class(sp_project),          intent(inout) :: self
        character(len=*), optional, intent(in)    :: dir
        character(len=:), allocatable :: ext, newname, stkname, abs_dir, nametmp, ext_out
        integer :: imic, nmics
        nmics = self%os_stk%get_noris()
        if( present(dir) )then
            call simple_mkdir( trim(dir), errmsg="sp_project::add_scale_tag" )
            abs_dir = simple_abspath( dir, 'sp_project :: add_scale_tag' )
        endif
        do imic=1,nmics
            call self%os_stk%getter(imic, 'stk', stkname)
            ext = fname2ext(trim(stkname))
            ext_out = '.'//ext
            if(present(dir))then
                nametmp = basename(add2fbody(stkname, '.'//trim(ext), trim(SCALE_SUFFIX)))
                newname = filepath(trim(abs_dir), trim(nametmp))
            else
                newname = add2fbody(stkname, '.'//trim(ext), trim(SCALE_SUFFIX))
            endif
            newname = fname_new_ext(newname, ext_out(2:))
            call self%os_stk%set(imic, 'stk', newname)
        end do
    end subroutine add_scale_tag

    ! os_out related methods

    subroutine add_cavgs2os_out( self, stk, smpd, imgkind )
        class(sp_project),          intent(inout) :: self
        character(len=*),           intent(in)    :: stk
        real,                       intent(in)    :: smpd ! sampling distance of images in stk
        character(len=*), optional, intent(in)    :: imgkind
        character(len=:), allocatable :: iimgkind
        character(len=LONGSTRLEN)     :: relpath
        integer                       :: ldim(3), nptcls, ind
        if( present(imgkind) )then
            allocate(iimgkind, source=trim(imgkind))
        else
            allocate(iimgkind, source='cavg')
        endif
        ! path and existence check
        call make_relativepath(CWD_GLOB, stk, relpath)
        ! find dimension of inputted stack
        call find_ldim_nptcls(relpath, ldim, nptcls)
        ! add os_out entry
        call self%add_entry2os_out(iimgkind, ind)
        ! fill-in field
        call self%os_out%set(ind, 'stk',     trim(relpath))
        call self%os_out%set(ind, 'box',     real(ldim(1)))
        call self%os_out%set(ind, 'nptcls',  real(nptcls))
        call self%os_out%set(ind, 'fromp',   1.0)
        call self%os_out%set(ind, 'top',     real(nptcls))
        call self%os_out%set(ind, 'smpd',    real(smpd))
        call self%os_out%set(ind, 'stkkind', 'single')
        call self%os_out%set(ind, 'imgkind', iimgkind)
        call self%os_out%set(ind, 'ctf',     'no')
        ! add congruent os_cls2D & os_cls3D
        if( self%os_cls2D%get_noris() /= nptcls )then
            call self%os_cls2D%new(nptcls, is_ptcl=.false.)
            call self%os_cls2D%set_all2single('state',1.)
        endif
        if( self%os_cls3D%get_noris() /= nptcls ) self%os_cls3D = self%os_cls2D
    end subroutine add_cavgs2os_out

    subroutine add_frcs2os_out( self, frc, which_imgkind )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: frc, which_imgkind
        character(len=LONGSTRLEN)        :: relpath
        integer                          :: ind
        select case(trim(which_imgkind))
            case('frc2D','frc3D')
                ! all good
            case DEFAULT
                THROW_HARD('invalid FRC kind: '//trim(which_imgkind)//'; add_frcs2os_out')
        end select
        ! full path and existence check
        call make_relativepath(CWD_GLOB, frc, relpath)
        ! add os_out entry
        call self%add_entry2os_out(which_imgkind, ind)
        ! fill-in field
        call self%os_out%set(ind, 'frcs',    trim(relpath))
        call self%os_out%set(ind, 'imgkind', trim(which_imgkind))
    end subroutine add_frcs2os_out

    subroutine add_fsc2os_out( self, fsc, state, box)
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fsc
        integer,           intent(in)    :: state, box
        character(len=:), allocatable :: imgkind
        character(len=LONGSTRLEN)     :: relpath
        integer                       :: i, ind, n_os_out
        ! full path and existence check
        call make_relativepath(CWD_GLOB, fsc, relpath)
        ! add os_out entry
        ! check if field is empty
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 )then
            n_os_out = 1
            ind      = 1
            call self%os_out%new(n_os_out, is_ptcl=.false.)
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
        call self%os_out%set(ind, 'fsc',     trim(relpath))
        call self%os_out%set(ind, 'imgkind', 'fsc')
        call self%os_out%set(ind, 'state',   real(state))
        call self%os_out%set(ind, 'box',     real(box))
    end subroutine add_fsc2os_out

    subroutine add_vol2os_out( self, vol, smpd, state, which_imgkind, box )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: vol, which_imgkind
        real,              intent(in)    :: smpd
        integer,           intent(in)    :: state
        integer, optional, intent(in)    :: box
        character(len=:), allocatable :: imgkind
        character(len=LONGSTRLEN)     :: relpath
        integer                       :: n_os_out, ind, i, ldim(3), ifoo
        select case(trim(which_imgkind))
            case('vol_cavg','vol','vol_msk')
                ! find_dimension of inputted volume
                call find_ldim_nptcls(vol, ldim, ifoo)
                if(present(box))then
                    if( ldim(1) /= box )then
                        THROW_HARD('invalid dimensions for volume: '//trim(vol)//'; add_vol2os_out 1')
                    endif
                endif
            case DEFAULT
                THROW_HARD('invalid VOL kind: '//trim(which_imgkind)//'; add_vol2os_out 3')
        end select
        ! path and existence check
        call make_relativepath(CWD_GLOB, vol, relpath)
        ! check if field is empty
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 )then
            n_os_out = 1
            ind      = 1
            call self%os_out%new(n_os_out, is_ptcl=.false.)
        else
            select case(trim(which_imgkind))
                case('vol_msk')
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
        call self%os_out%set(ind, 'vol',     trim(relpath))
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
            call self%os_out%new(n_os_out, is_ptcl=.false.)
        else
            ind = 0
            do i=1,n_os_out
                if( self%os_out%isthere(i,'imgkind') )then
                    imgkind = trim(self%os_out%get_static(i,'imgkind'))
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

    subroutine get_cavgs_stk( self, stkname, ncls, smpd, imgkind, fail )
        class(sp_project),             intent(inout) :: self
        character(len=:), allocatable, intent(inout) :: stkname
        integer,                       intent(out)   :: ncls
        real,                          intent(out)   :: smpd
        character(len=*), optional,    intent(in)    :: imgkind
        logical,          optional,    intent(in)    :: fail
        character(len=:), allocatable :: ikind, iimgkind
        integer :: n_os_out, ind, i, cnt
        logical :: fail_here
        if( present(imgkind) )then
            allocate(iimgkind, source=trim(imgkind))
        else
            allocate(iimgkind, source='cavg')
        endif
        fail_here = .true.
        if( present(fail) )fail_here = fail
        ! check if field is empty
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 )then
            if( fail_here )then
                THROW_HARD('trying to fetch from empty os_out field; get_cavgs_stk')
            else
                stkname = NIL
                ncls    = 0
                smpd    = 0.
                return
            endif
        endif
        ! look for cavgs
        ind = 0
        cnt = 0
        do i=1,n_os_out
            if( self%os_out%isthere(i,'imgkind') )then
                ikind = trim(self%os_out%get_static(i,'imgkind'))
                if(trim(ikind).eq.trim(iimgkind))then
                    ind = i
                    cnt = cnt + 1
                endif
            endif
        end do
        if( fail_here )then
            if( cnt > 1 )  THROW_HARD('multiple os_out entries with imgkind='//iimgkind//', aborting... get_cavgs_stk')
            if( cnt == 0 ) THROW_HARD('no os_out entry with imgkind='//iimgkind//' identified, aborting... get_cavgs_stk')
        else if( cnt > 1 .or. cnt == 0 )then
            stkname = NIL
            ncls    = 0
            smpd    = 0.
            return
        endif
        ! set return values
        if( allocated(stkname) ) deallocate(stkname)
        stkname = trim(self%os_out%get_static(ind,'stk'))
        ncls = nint(self%os_out%get(ind, 'nptcls'))
        smpd = self%os_out%get(ind, 'smpd')
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
            case('vol_cavg','vol','vol_msk')
                ! all good
            case DEFAULT
                THROW_HARD('invalid VOL kind: '//trim(imgkind)//'; get_vol')
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
        case('vol_msk')
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
                THROW_HARD('no os_out entry with imgkind=volXXX identified, aborting...; get_vol')
            endif
        endif
        if( cnt > 1 )  THROW_HARD('multiple os_out entries with imgkind=volXXX, aborting...; get_vol')
        ! set output
        deallocate(vol_fname)
        call self%os_out%getter(ind, 'vol', vol_fname)
        smpd = self%os_out%get(ind,  'smpd')
        box  = nint(self%os_out%get(ind,  'box'))
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
        if( cnt == 0 )THROW_HARD('no os_out entry with imgkind=fsc identified, aborting...; get_fsc')
        if( cnt > 1 ) THROW_HARD('multiple os_out entries with imgkind=fsc, aborting...; get_fsc')
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
        logical :: fail_here, found
        select case(trim(which_imgkind))
            case('frc2D','frc3D')
                ! all good
            case DEFAULT
                THROW_HARD('invalid FRC kind: '//trim(which_imgkind)//'; get_frcs')
        end select
        fail_here = .true.
        if( present(fail) )fail_here = fail
        ! check if field is empty
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 ) THROW_HARD('trying to fetch from empty os_out field; get_frcs')
        ! look for cavgs
        ind   = 0
        cnt   = 0
        found = .false.
        do i=1,n_os_out
            if( self%os_out%isthere(i,'imgkind') )then
                call self%os_out%getter(i,'imgkind',imgkind)
                if(trim(imgkind).eq.trim(which_imgkind))then
                    ind   = i
                    cnt   = cnt + 1
                    found = .true.
                endif
            endif
        end do
        if( allocated(frcs) ) deallocate(frcs)
        if( fail_here )then
            if( cnt > 1 )  THROW_HARD('multiple os_out entries with imgkind=frcXD, aborting...; get_frcs')
            if( cnt == 0 ) THROW_HARD('no os_out entry with imgkind=frcsXD identified, aborting...; get_frcs')
        endif
        ! set return values
        if( found )then
            call self%os_out%getter(ind,'frcs',frcs)
        else
            frcs = NIL
        endif
    end subroutine get_frcs

    subroutine get_imginfo_from_osout( self, smpd, box, nptcls )
        class(sp_project), intent(inout) :: self
        real,              intent(out)   :: smpd
        integer,           intent(out)   :: box, nptcls
        character(len=LONGSTRLEN) :: imgkind
        integer :: i, nos
        nptcls = 0
        smpd   = 0.
        box    = 0
        nos    = self%os_out%get_noris()
        if( nos == 0 )return
        ! nptcls: look for cavgs, defaults to zero for other entries (frcs, volumes)
        do i=1,nos
            if( self%os_out%isthere(i,'imgkind') )then
                imgkind = trim(self%os_out%get_static(i,'imgkind'))
                if(trim(imgkind).eq.'cavg')then
                    if( self%os_out%isthere(i,'fromp').and.self%os_out%isthere(i,'top') )then
                        nptcls = nptcls + nint(self%os_out%get(i,'top')) - nint(self%os_out%get(i,'fromp')) + 1
                    else
                        THROW_HARD('Missing fromp and top entries in cavg; get_imginfo_from_osout')
                    endif
                endif
            endif
        end do
        ! box/smpd: first in
        do i=1,nos
            if( self%os_out%isthere(i,'smpd').and.self%os_out%isthere(i,'box') )then
                smpd = self%os_out%get(i,'smpd')
                box  = nint(self%os_out%get(i,'box'))
                return
            endif
        end do
    end subroutine get_imginfo_from_osout

    ! getters

    integer function get_n_insegment( self, oritype )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        class(oris), pointer :: pos => NULL()
        get_n_insegment = 0
        select case(oritype)
            case('stk')
                get_n_insegment = self%get_nstks()
            case('ptcl2D','ptcl3D')
                ! # ptcl2D = # ptcl3D
                get_n_insegment = self%get_nptcls()
                if( get_n_insegment /= self%os_ptcl2D%get_noris() .or.&
                   &get_n_insegment /= self%os_ptcl3D%get_noris() )then
                   THROW_HARD('Inconstitent number of particles in STK/PTCL2D/PTCL3D segments; get_n_insegment')
                endif
            case DEFAULT
                call self%ptr2oritype(oritype, pos)
                get_n_insegment = pos%get_noris()
                nullify(pos)
        end select
    end function get_n_insegment

    integer function get_n_insegment_state( self, oritype, state )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        real,                      intent(in)    :: state
        class(oris), pointer :: pos => NULL()
        integer              :: iori
        get_n_insegment_state = 0
        select case(oritype)
            case('ptcl2D','ptcl3D')
                ! # ptcl2D = # ptcl3D
                if( self%os_ptcl2D%get_noris() /= self%os_ptcl3D%get_noris() )then
                   THROW_HARD('Inconstitent number of particles in STK/PTCL2D/PTCL3D segments; get_n_insegment_state')
                endif
        end select
        call self%ptr2oritype(oritype, pos)
        do iori=1,pos%get_noris()
            if(.not. pos%isthere(iori,'state') )then
                THROW_HARD('state flag missing from ori; get_n_insegment_state')
            endif
            if( pos%get_state(iori) == state )then
                get_n_insegment_state = get_n_insegment_state + 1
            endif
        enddo
        nullify(pos)
    end function get_n_insegment_state

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
                write(logfhandle,*) 'nptcls from ptcls', get_nptcls
                write(logfhandle,*) 'nptcls from top', nint(self%os_stk%get(nos,'top'))
                THROW_HARD('total # particles .ne. last top index; get_nptcls')
            endif
        endif
    end function get_nptcls

    integer function get_box( self )
        class(sp_project), target, intent(inout) :: self
        integer :: n_os_stk
        get_box  = 0
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk == 0 )then
            THROW_HARD('empty os_stk field! get_box')
        endif
        get_box = nint( self%os_stk%get(1,'box') )
    end function get_box

    logical function has_boxcoords(self, iptcl)
        class(sp_project), target, intent(in) :: self
        integer,                   intent(in) :: iptcl
        has_boxcoords = .false.
        if( self%os_ptcl2D%isthere(iptcl,'xpos') .and. self%os_ptcl2D%isthere(iptcl,'ypos'))then
            has_boxcoords = .true.
        endif
    end function has_boxcoords

    subroutine get_boxcoords( self, iptcl, coords )
        class(sp_project), target, intent(in)  :: self
        integer,                   intent(in)  :: iptcl
        integer,                   intent(out) :: coords(2)
        coords = 0
        if( self%has_boxcoords(iptcl) )then
            coords(1) = nint(self%os_ptcl2D%get(iptcl,'xpos'))
            coords(2) = nint(self%os_ptcl2D%get(iptcl,'ypos'))
        endif
    end subroutine get_boxcoords

    real function get_smpd( self )
        class(sp_project), target, intent(inout) :: self
        integer :: n_os_stk
        get_smpd  = 0.
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk == 0 )then
            THROW_HARD('empty os_stk field! get_smpd')
        endif
        get_smpd = self%os_stk%get(1,'smpd')
    end function get_smpd

    integer function get_nmovies( self )
        class(sp_project), target, intent(inout) :: self
        character(len=:), allocatable :: imgkind
        integer :: i
        get_nmovies = 0
        do i=1,self%os_mic%get_noris()
            call self%os_mic%getter(i,'imgkind',imgkind)
            if( trim(imgkind).eq.'movie' .or. trim(imgkind).eq.'mic' ) get_nmovies = get_nmovies + 1
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

    integer function get_nframes( self )
        class(sp_project), target, intent(inout) :: self
        integer :: i
        get_nframes = 0
        do i=1,self%os_mic%get_noris()
            if( self%os_mic%isthere(i,'frame') )get_nframes = get_nframes + 1
        enddo
    end function get_nframes

    integer function get_nstks( self )
        class(sp_project), target, intent(in) :: self
        get_nstks = self%os_stk%get_noris()
    end function get_nstks

    character(len=STDLEN) function get_ctfflag( self, oritype, iptcl )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,         optional, intent(in)    :: iptcl
        class(oris), pointer         :: ptcl_field
        character(len=:),allocatable :: ctfflag_str
        integer :: stkind, ind_in_stk, ind
        get_ctfflag = 'no'
        nullify(ptcl_field)
        ! set field pointer
        select case(trim(oritype))
            case('ptcl2D')
                ptcl_field => self%os_ptcl2D
            case('ptcl3D')
                ptcl_field => self%os_ptcl3D
            case('cls2D', 'cls3D')
                return
            case DEFAULT
                THROW_HARD('oritype: '//trim(oritype)//' not supported by get_ctfflag')
        end select
        ! do the index mapping
        ind = 1
        if( present(iptcl) ) ind = iptcl
        call self%map_ptcl_ind2stk_ind(oritype, ind, stkind, ind_in_stk)
        ! CTF flag string
        if( self%os_stk%isthere(stkind, 'ctf') )then
            call self%os_stk%getter(stkind, 'ctf', ctfflag_str)
        else if( ptcl_field%isthere(ind, 'ctf') )then
            call ptcl_field%getter(ind, 'ctf', ctfflag_str)
        else
           ctfflag_str = 'no'
        endif
        get_ctfflag = trim(ctfflag_str)
    end function get_ctfflag

    integer(kind(ENUM_CTFFLAG)) function get_ctfflag_type( self, oritype, iptcl )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,         optional, intent(in)    :: iptcl
        character(len=:), allocatable :: ctfflag
        ctfflag = self%get_ctfflag(oritype, iptcl=iptcl)
        select case(trim(ctfflag))
            case('no')
                get_ctfflag_type = CTFFLAG_NO
            case('yes')
                get_ctfflag_type = CTFFLAG_YES
            case('mul')
                THROW_HARD('ctf=mul depreciated; get_ctfflag_type')
            case('flip')
                get_ctfflag_type = CTFFLAG_FLIP
            case DEFAULT
                THROW_HARD('unsupported ctf flag: '//trim(ctfflag)//'; get_ctfflag_type')
        end select
    end function get_ctfflag_type

    ! the entire project must be phase-plate, so the 1 is
    logical function has_phaseplate( self, oritype, iptcl )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,         optional, intent(in)    :: iptcl
        character(len=:), allocatable :: phaseplate
        integer :: ind
        has_phaseplate = .false.
        if( trim(oritype) .eq. 'cls3D' ) return
        ind = 1
        if( present(iptcl) ) ind = iptcl
        ! get info
        if( self%os_stk%isthere(ind, 'phaseplate') )then
            phaseplate = trim(self%os_stk%get_static(ind, 'phaseplate'))
        else
            phaseplate = 'no'
        endif
        has_phaseplate = trim(phaseplate).eq.'yes'
    end function has_phaseplate

    logical function has_boxfile( self )
        class(sp_project), target, intent(inout) :: self
        has_boxfile = self%os_mic%isthere('boxfile')
    end function has_boxfile

    function get_ctfparams( self, oritype, iptcl ) result( ctfvars )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,                   intent(in)    :: iptcl
        class(oris), pointer  :: ptcl_field
        character(len=STDLEN) :: ctfflag, phaseplate
        type(ctfparams)       :: ctfvars
        real                  :: smpd
        integer               :: stkind, ind_in_stk, box, ncls
        logical               :: dfy_was_there, l_noctf
        nullify(ptcl_field)
        ! set field pointer
        select case(trim(oritype))
            case('stk')
                ptcl_field => self%os_stk
                stkind = iptcl
            case('ptcl2D')
                ptcl_field => self%os_ptcl2D
                call self%map_ptcl_ind2stk_ind(oritype, iptcl, stkind, ind_in_stk)
            case('ptcl3D')
                ptcl_field => self%os_ptcl3D
                call self%map_ptcl_ind2stk_ind(oritype, iptcl, stkind, ind_in_stk)
            case('cls3D')
                call self%get_imginfo_from_osout(smpd, box, ncls)
                ctfvars%ctfflag = CTFFLAG_NO
                ctfvars%smpd    = smpd
                ctfvars%kv      = 0.
                ctfvars%cs      = 0.
                ctfvars%fraca   = 0.
                ctfvars%dfx     = 0.
                ctfvars%dfy     = 0.
                ctfvars%angast  = 0.
                ctfvars%phshift = 0.
                ctfvars%l_phaseplate = .false.
                return
            case DEFAULT
                THROW_HARD('oritype: '//trim(oritype)//' not supported by get_ctfparams')
        end select
        ! extract the CTF parameters
        ! sampling distance
        if( self%os_stk%isthere(stkind, 'smpd') )then
            ctfvars%smpd = self%os_stk%get(stkind, 'smpd')
        else
            THROW_HARD('smpd (sampling distance) lacking in os_stk and ptcl fields; get_ctfparams')
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
                THROW_HARD('ctf key lacking in os_stk_field & ptcl fields; get_ctfparams')
            case('no')
                ctfvars%ctfflag = CTFFLAG_NO
                l_noctf = .true.
            case('yes')
                ctfvars%ctfflag = CTFFLAG_YES
            case('mul')
                THROW_HARD('ctf=mul depreciated; get_ctfparams')
            case('flip')
                ctfvars%ctfflag = CTFFLAG_FLIP
            case DEFAULT
                write(logfhandle,*) 'stkind/iptcl: ', stkind, iptcl
                THROW_HARD('unsupported ctf flag: '// trim(ctfflag)//'; get_ctfparams')
        end select
        ! acceleration voltage
        if( self%os_stk%isthere(stkind, 'kv') )then
            ctfvars%kv = self%os_stk%get(stkind, 'kv')
        else
            THROW_HARD('kv (acceleration voltage) lacking in os_stk_field; get_ctfparams')
        endif
        ! spherical aberration constant
        if( self%os_stk%isthere(stkind, 'cs') )then
            ctfvars%cs = self%os_stk%get(stkind, 'cs')
        else
            THROW_HARD('cs (spherical aberration constant) lacking in os_stk_field; get_ctfparams')
        endif
        ! fraction of amplitude contrast
        if( self%os_stk%isthere(stkind, 'fraca') )then
            ctfvars%fraca = self%os_stk%get(stkind, 'fraca')
        else
            THROW_HARD('fraca (fraction of amplitude contrast) lacking in os_stk_field; get_ctfparams')
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
            ctfvars%dfx = ptcl_field%get_dfx(iptcl)
        else
            call ptcl_field%print_(iptcl)
            THROW_HARD('dfx (defocus in x) lacking in ptcl_field; get_ctfparams')
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
                print *, 'iptcl: ', iptcl
                print *, 'oritype in get_ctfparams ', trim(oritype)
                call ptcl_field%print_(iptcl)
                THROW_HARD('astigmatic CTF model requires angast (angle of astigmatism) lacking in os_stk field; get_ctfparams')
            else
                ctfvars%angast = 0.
            endif
        endif
        ! has phaseplate
        if( self%os_stk%isthere(stkind, 'phaseplate') )then
            phaseplate = trim(self%os_stk%get_static(stkind, 'phaseplate'))
        else
            phaseplate = 'no'
        endif
        ctfvars%l_phaseplate = trim(phaseplate).eq.'yes'
        ! additional phase shift
        ctfvars%phshift = 0.
        if( ptcl_field%isthere(iptcl, 'phshift') )then
            ctfvars%phshift = ptcl_field%get(iptcl, 'phshift')
        endif
    end function get_ctfparams

    subroutine get_sp_oris( self, which_imgkind, os )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which_imgkind
        class(oris),       intent(inout) :: os
        select case(trim(which_imgkind))
            case('mic')
                call os%copy(self%os_mic,    is_ptcl=.false.)
            case('stk')
                call os%copy(self%os_stk,    is_ptcl=.false.)
            case('ptcl2D')
                call os%copy(self%os_ptcl2D, is_ptcl=.true.)
            case('cls2D')
                call os%copy(self%os_cls2D,  is_ptcl=.false.)
            case('cls3D')
                call os%copy(self%os_cls3D,  is_ptcl=.false.)
            case('ptcl3D')
                call os%copy(self%os_ptcl3D, is_ptcl=.true.)
            case('out')
                call os%copy(self%os_out,    is_ptcl=.false.)
            case('projinfo')
                call os%copy(self%projinfo,  is_ptcl=.false.)
            case('jobproc')
                call os%copy(self%jobproc,   is_ptcl=.false.)
            case('compenv')
                call os%copy(self%compenv,   is_ptcl=.false.)
            case DEFAULT
                THROW_HARD('unsupported which_imgkind flag; get_sp_oris')
        end select
    end subroutine get_sp_oris

    subroutine ptr2oritype( self, oritype, os_ptr )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        class(oris),      pointer, intent(inout) :: os_ptr
        select case(trim(oritype))
            case('mic')
                os_ptr => self%os_mic
            case('stk')
                os_ptr => self%os_stk
            case('ptcl2D')
                os_ptr => self%os_ptcl2D
            case('cls2D')
                os_ptr => self%os_cls2D
            case('cls3D')
                os_ptr => self%os_cls3D
            case('ptcl3D')
                os_ptr => self%os_ptcl3D
            case('out')
                os_ptr => self%os_out
            case DEFAULT
                THROW_HARD('unsupported oritype: '//trim(oritype)//'; ptr2segment')
        end select
    end subroutine ptr2oritype

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
                THROW_HARD('oritype: '//trim(oritype)//' not supported by is_virgin_field')
        end select
        n = os%get_noris()
        if( n == 0 )then
            THROW_WARN('cannot check virginity of non-existent field (touched for the very first time???)')
            return
        endif
        do i=1,n
            if( os%has_been_searched(i) )return
        enddo
        is_virgin_field = .true.
    end function is_virgin_field

    subroutine get_mic2stk_inds( self, mic2stk_inds, stk2mic_inds )
        class(sp_project),    intent(inout) :: self
        integer, allocatable, intent(inout) :: mic2stk_inds(:), stk2mic_inds(:)
        integer :: imic,istk,nmics,nstks,nptcls_mic,nptcls_stk,state_mic,state_stk
        if(allocated(mic2stk_inds))deallocate(mic2stk_inds)
        if(allocated(stk2mic_inds))deallocate(stk2mic_inds)
        nmics = self%os_mic%get_noris()
        nstks = self%os_stk%get_noris()
        if( nmics==0 .or. nstks==0 )then
            THROW_WARN('Empty fields! Fields need be populated; get_mic2stk_inds')
        endif
        if( nmics < nstks )then
            THROW_HARD('MIC & STK fileds indexing error 1! get_mic2stk_inds')
        endif
        allocate(mic2stk_inds(nmics),stk2mic_inds(nstks),source=0)
        if( nmics == nstks )then
            do imic = 1,nmics
                mic2stk_inds(imic) = imic
                stk2mic_inds(imic) = imic
            enddo
        else
            istk = 0
            do imic = 1,nmics
                nptcls_mic = nint(self%os_mic%get(imic,'nptcls'))
                if( nptcls_mic > 0 )then
                    istk = istk+1
                    if( istk > nstks ) THROW_HARD('Too many stacks!  get_mic2stk_inds')
                else
                    ! micrographs without particles have no stack
                    cycle
                endif
                mic2stk_inds(imic) = istk
                stk2mic_inds(istk) = imic
            enddo
        endif
        ! state consistency
        do imic = 1,nmics
            istk = mic2stk_inds(imic)
            if( istk == 0 ) cycle
            nptcls_mic = nint(self%os_mic%get(imic,'nptcls'))
            nptcls_stk = nint(self%os_stk%get(istk,'nptcls'))
            if( nptcls_mic /= nptcls_stk )then
                 THROW_HARD('Inconsistent number of particles!  get_mic2stk_inds')
            endif
            state_mic = self%os_mic%get_state(imic)
            state_stk = self%os_stk%get_state(istk)
            if( state_mic /= state_stk )then
                THROW_HARD('Inconsistent state!  get_mic2stk_inds')
            endif
        enddo
    end subroutine get_mic2stk_inds

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
                THROW_HARD('unsupported which_imgkind flag; set_sp_oris')
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
        real    :: scale_factor, smpd_sc, smpd
        integer :: box, box_sc, istk, n_os_stk
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk == 0 ) THROW_HARD('Empty stack object! scale_projfile')
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
        if( cline%defined('mkdir') )     call cline_scale%set('mkdir', cline%get_carg('mkdir'))
        if( present(dir) )               call cline_scale%set('dir_target',trim(dir)//path_separator)
        if( box == box_sc )then
            ! no scaling
            new_projfile = trim(projfile)
            return
        endif
        ! parameter updates
        do istk = 1,n_os_stk
            call self%os_stk%set(istk, 'smpd', real(smpd_sc))
            call self%os_stk%set(istk, 'box', real(box_sc))
        enddo
        call self%os_ptcl2D%mul_shifts(scale_factor)
        call self%os_ptcl3D%mul_shifts(scale_factor)
        ! name changes and list for scaling job
        new_projname = trim(projname)//SCALE_SUFFIX
        new_projfile = trim(new_projname)//trim(METADATA_EXT)
        call cline%set('projname', trim(new_projname))
        call cline%delete('projfile')
        call self%update_projinfo(cline)
        if(present(dir))then
            call self%add_scale_tag(dir=trim(dir)//path_separator)
        else
            call self%add_scale_tag
        endif
        ! save
        call self%write
        ! command line for scaling
        call cline_scale%set('newbox', real(box_sc))
        if( cline%defined('state') )  call cline_scale%set('state', cline%get_rarg('state'))
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
        integer,          allocatable :: parts(:,:)
        character(len=:), allocatable :: fname, projfile
        type(str4arr),    allocatable :: os_strings(:)
        class(oris),          pointer :: os => null()
        type(binoris) :: bos_doc
        integer       :: i, numlen, n_records, partsz, isegment, strlen_max
        numlen = len(int2str(ndocs))
        if( present(numlen_in) ) numlen = numlen_in
        parts  = split_nobjs_even(nptcls, ndocs)
        ! convert from flag to enumerator
        isegment = oritype2segment(oritype)
        if( isegment==3 .or. isegment==6 )then
            call self%ptr2oritype(oritype, os)
            do i=1,ndocs
                ! read part
                fname     = trim(adjustl(fbody))//int2str_pad(i,numlen)//'.simple'
                call bos_doc%open(trim(fname))
                n_records = bos_doc%get_n_records(isegment)
                partsz    = parts(i,2) - parts(i,1) + 1
                if( n_records /= partsz )then
                    write(logfhandle,*) 'ERROR, # records does not match expectation'
                    write(logfhandle,*) 'EXTRACTED FROM file: ', trim(fname)
                    write(logfhandle,*) 'n_records: ', n_records
                    write(logfhandle,*) 'CALCULATED FROM input p%nptcls/p%ndocs'
                    write(logfhandle,*) 'fromto: ', parts(i,1), parts(i,2)
                    write(logfhandle,*) 'partsz: ', partsz
                    stop
                endif
                call bos_doc%read_segment(isegment, os)
                call bos_doc%close
            end do
            ! write
            call self%projinfo%getter(1, 'projfile', projfile)
            call self%bos%open(projfile)
            call self%bos%write_segment_inside(isegment, os)
        else
            ! allocate merged string representation
            allocate(os_strings(nptcls))
            ! maxium string length
            strlen_max = 0
            ! read into string representation
            do i=1,ndocs
                ! read part
                fname     = trim(adjustl(fbody))//int2str_pad(i,numlen)//'.simple'
                call bos_doc%open(trim(fname))
                n_records = bos_doc%get_n_records(isegment)
                partsz    = parts(i,2) - parts(i,1) + 1
                if( n_records /= partsz )then
                    write(logfhandle,*) 'ERROR, # records does not match expectation'
                    write(logfhandle,*) 'EXTRACTED FROM file: ', trim(fname)
                    write(logfhandle,*) 'n_records: ', n_records
                    write(logfhandle,*) 'CALCULATED FROM input p%nptcls/p%ndocs'
                    write(logfhandle,*) 'fromto: ', parts(i,1), parts(i,2)
                    write(logfhandle,*) 'partsz: ', partsz
                    stop
                endif
                call bos_doc%read_segment(isegment, os_strings)
                strlen_max = max(strlen_max, bos_doc%get_n_bytes_per_record(isegment))
                call bos_doc%close
            end do
            ! write
            call self%projinfo%getter(1, 'projfile', projfile)
            call self%bos%open(projfile)
            call self%bos%write_segment_inside(isegment, os_strings, [1,nptcls], strlen_max)
            ! transfer to memory & destruct
            call self%ptr2oritype(oritype, os)
            do i=1,nptcls
                call os%str2ori(i, os_strings(i)%str)
                if( allocated(os_strings(i)%str) ) deallocate(os_strings(i)%str)
            end do
            deallocate(os_strings)
        endif
        nullify(os)
        ! no need to update header (taken care of in binoris object)
        call self%bos%close
    end subroutine merge_algndocs

    ! map shifts obtained by cluster2D to the 3D field for cases where an external
    ! starting model is used for initializing 3D refinement (nanoparticles)
    subroutine map2Dshifts23D( self )
        class(sp_project), intent(inout) :: self
        integer :: noris_ptcl3D, noris_ptcl2D
        real, allocatable :: shifts(:)
        if( self%is_virgin_field('ptcl2D') )then
            return
        endif
        ! ensure ptcl3D field congruent with ptcl2D field
        noris_ptcl3D = self%os_ptcl3D%get_noris()
        noris_ptcl2D = self%os_ptcl2D%get_noris()
        if( noris_ptcl3D /= noris_ptcl2D )then
            ! preserve defocus parameters, stack indices
            self%os_ptcl3D = self%os_ptcl2D
            call self%os_ptcl3D%delete_2Dclustering(keepshifts=.true., keepcls=.true.)
        else
            ! transfer shifts
            shifts = self%os_ptcl2D%get_all('x')
            call self%os_ptcl3D%set_all('x', shifts)
            deallocate(shifts)
            shifts = self%os_ptcl2D%get_all('y')
            call self%os_ptcl3D%set_all('y', shifts)
            deallocate(shifts)
        endif
    end subroutine map2Dshifts23D

    ! this map2ptcls routine assumes that any selection of class averages is done
    ! exclusively by state=0 flagging without any physical deletion
    subroutine map2ptcls( self )
        class(sp_project), intent(inout) :: self
        integer, allocatable :: particles(:)
        type(ori) :: ori2d, ori_comp, o
        integer   :: ncls, icls, iptcl, pind, noris_ptcl3D, noris_ptcl2D
        real      :: corr, rproj, rstate
        if( self%is_virgin_field('cls3D') )then
            THROW_HARD('os_cls3D is virgin field; nothing to map back; map2ptcls')
        endif
        ! in case 2D was not done with SIMPLE but class averages imported from elsewhere
        if( self%os_ptcl2D%get_noris() == 0 ) return
        if( self%is_virgin_field('ptcl2D')  ) return
        ! ensure ptcl3D field congruent with ptcl2D field
        noris_ptcl3D = self%os_ptcl3D%get_noris()
        noris_ptcl2D = self%os_ptcl2D%get_noris()
        if( noris_ptcl3D /= noris_ptcl2D )then
            ! preserve defocus parameters, stack indices
            self%os_ptcl3D = self%os_ptcl2D
            call self%os_ptcl3D%delete_2Dclustering(keepcls=.true.)
        endif
        ! do the mapping
        ncls = self%os_cls3D%get_noris()
        do icls=1,ncls
            ! get particle indices
            call self%os_ptcl2D%get_pinds(icls, 'class', particles)
            if( allocated(particles) )then
                ! get 3d ori info
                call self%os_cls3D%get_ori(icls, o)
                rproj  = o%get('proj')
                rstate = o%get('state')
                corr   = o%get('corr')
                do iptcl=1,size(particles)
                    ! get particle index
                    pind  = particles(iptcl)
                    ! get 2d ori
                    call self%os_ptcl2D%get_ori(pind, ori2d)
                    ! transfer original parameters in self%os_ptcl2D
                    call self%os_ptcl2D%get_ori(pind, ori_comp)
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
        call ori2d%kill
        call ori_comp%kill
        call o%kill
    end subroutine map2ptcls

    ! this map2ptcls routine assumes that any selection of class averages is done
    ! exclusively by state=0 flagging without any physical deletion
    subroutine map2ptcls_state( self )
        class(sp_project), intent(inout) :: self
        integer, allocatable :: particles(:)
        integer   :: ncls, icls, iptcl, pind, noris_ptcl3D, noris_ptcl2D
        real      :: rstate
        noris_ptcl2D = self%os_ptcl2D%get_noris()
        if( noris_ptcl2D == 0 )then
            THROW_WARN('empty PTCL2D field. Nothing to do; map2ptcls_state')
            return
        endif
        ! ensure ptcl3D field congruent with ptcl2D field
        noris_ptcl3D = self%os_ptcl3D%get_noris()
        if( noris_ptcl3D /= noris_ptcl2D )then
            ! preserve defocus parameters, stack indices
            self%os_ptcl3D = self%os_ptcl2D
            call self%os_ptcl3D%delete_2Dclustering(keepcls=.true.)
        endif
        ! undo previous selection & excludes non classified particles
        do iptcl=1,noris_ptcl2D
            if( .not.self%os_ptcl2D%isthere(iptcl, 'class') )then
                call self%os_ptcl2D%set(iptcl, 'state', 0.)
                call self%os_ptcl3D%set(iptcl, 'state', 0.)
            else
                 call self%os_ptcl2D%set(iptcl, 'state', 1.)
                 call self%os_ptcl3D%set(iptcl, 'state', 1.)
            endif
        end do
        ! do the mapping
        ncls = self%os_cls2D%get_noris()
        do icls=1,ncls
            ! get particle indices
            call self%os_ptcl2D%get_pinds(icls, 'class', particles)
            if( allocated(particles) )then
                ! get 3d ori info
                rstate = self%os_cls2D%get(icls,'state')
                do iptcl=1,size(particles)
                    ! get particle index
                    pind = particles(iptcl)
                    ! update state in self%os_ptcl2D
                    call self%os_ptcl2D%set(pind, 'state', rstate)
                    ! update state in self%os_ptcl3D
                    call self%os_ptcl3D%set(pind, 'state', rstate)
                end do
                deallocate(particles)
            endif
        end do
        ! cls3D mirrors cls2D
        if( self%os_cls3D%get_noris() == ncls)then
            do icls=1,ncls
                call self%os_cls3D%set(icls, 'state', self%os_cls2D%get(icls,'state'))
            enddo
        else if( self%os_cls3D%get_noris() > 0 )then
            THROW_WARN('Inconsistent number of classes in cls2D & cls3D segments')
        endif
        ! state = 0 all entries that don't have a state/class label
        do iptcl=1,noris_ptcl2D
            if( .not.self%os_ptcl2D%isthere(iptcl, 'state') .or.&
                &.not.self%os_ptcl3D%isthere(iptcl, 'state') )then
                 call self%os_ptcl2D%set(iptcl, 'state', 0.)
                 call self%os_ptcl3D%set(iptcl, 'state', 0.)
            endif
        end do
    end subroutine map2ptcls_state

    subroutine replace_project( self, projfile_src, oritype )
        class(sp_project),     intent(inout) :: self
        character(len=*),      intent(in)    :: projfile_src
        character(len=STDLEN), intent(in)    :: oritype
        character(len=:), allocatable :: boxfname, stkfname, movfname, micfname, src_path
        character(len=LONGSTRLEN)     :: relstkname, relboxname, relmicname, relmovname
        type(sp_project) :: self_src
        type(ori)        :: o, o_src
        integer          :: istk, nstks, imic, nmics, iptcl
        logical          :: err
        src_path = get_fpath(projfile_src)
        call self_src%read(projfile_src)
        select case(trim(oritype))
        case('ptcl2D','ptcl3D')
            if(  self%os_ptcl2D%get_noris() /= self_src%os_ptcl2D%get_noris() .or.&
                &self%os_ptcl3D%get_noris() /= self_src%os_ptcl3D%get_noris())then
                THROW_HARD('Inconsistent # of ptcls in project files!')
            endif
            if( trim(oritype) == 'ptcl2D' )then
                if( self_src%os_ptcl2D%get_noris() == 0 ) return
                ! transfer eo, state, weight and alignement parameters. Defocus untouched
                do iptcl=1,self%os_ptcl2D%get_noris()
                    call self_src%os_ptcl2D%get_ori(iptcl, o_src)
                    call self%os_ptcl2D%transfer_2Dparams(iptcl, o_src)
                    call self%os_ptcl2D%set(iptcl,'state', self_src%os_ptcl2D%get(iptcl,'state'))
                    call self%os_ptcl3D%set(iptcl,'state', self_src%os_ptcl2D%get(iptcl,'state'))
                enddo
            else
                if( self_src%os_ptcl3D%get_noris() == 0 ) return
                ! transfer eo, state, weight and alignement parameters. Defocus untouched
                do iptcl=1,self%os_ptcl3D%get_noris()
                    call self_src%os_ptcl3D%get_ori(iptcl, o_src)
                    call self%os_ptcl3D%transfer_3Dparams(iptcl, o_src)
                    call self%os_ptcl3D%set(iptcl,'state', self_src%os_ptcl3D%get(iptcl,'state'))
                    call self%os_ptcl2D%set(iptcl,'state', self_src%os_ptcl3D%get(iptcl,'state'))
                enddo
            endif
        case('stk')
            if( self%os_stk%get_noris() /= self_src%os_stk%get_noris())then
                THROW_HARD('Inconsistent # of stk/ptcls in project files!')
            endif
            nstks = self_src%os_stk%get_noris()
            if( nint(self_src%os_stk%get(nstks,'top')) /= self%os_ptcl2D%get_noris() )then
                THROW_HARD('Inconsistent # of particles between stk & particles fields')
            endif
            if( nstks == 0 ) return
            do istk = 1,nstks
                err   = .false.
                call self_src%os_stk%get_ori(istk, o_src)
                call self%os_stk%get_ori(istk, o)
                if( nint(o%get('box')) /= nint(o_src%get('box')) )then
                    THROW_HARD('Inconsistent box size')
                endif
                if( .not.is_equal(o%get('smpd'),o_src%get('smpd')) )then
                    THROW_HARD('Inconsistent sampling distance')
                endif
                if( o_src%isthere('stk') )then
                    call o_src%getter('stk',stkfname)
                    if( stkfname(1:1).ne.'/' ) stkfname = trim(src_path)//'/'//trim(stkfname)
                    call make_relativepath(CWD_GLOB, stkfname, relstkname)
                    if( file_exists(relstkname) )then
                        call o_src%set('stk', relstkname)
                    else
                        err = .true.
                    endif
                else
                    err = .true.
                endif
                if( o_src%isthere('boxfile') )then
                    call o_src%getter('boxfile',boxfname)
                    if( boxfname(1:1).ne.'/' ) boxfname = trim(src_path)//'/'//trim(boxfname)
                    call make_relativepath(CWD_GLOB, boxfname, relboxname)
                    if( file_exists(relboxname) )then
                        call o_src%set('boxfile', relboxname)
                    else
                        err = .true.
                    endif
                endif
                if( err )then
                    call o_src%print_ori
                    write(logfhandle,*) trim(relstkname)
                    write(logfhandle,*) trim(relboxname)
                    THROW_HARD('Missing stack or boxfile!')
                else
                    call self%os_stk%set_ori(istk,o_src)
                endif
            enddo
        case('mic')
            if( self%os_mic%get_noris() /= self_src%os_mic%get_noris())then
                THROW_HARD('Inconsistent # of movies/micrographs in project files!')
            endif
            nmics = self_src%os_mic%get_noris()
            if( nmics == 0 ) return
            do imic = 1,nmics
                call self%os_mic%get_ori(imic, o)
                call self_src%os_mic%get_ori(imic, o_src)
                if(  nint(o%get('xdim')) /= nint(o_src%get('xdim')) .or.&
                    &nint(o%get('ydim')) /= nint(o_src%get('ydim')))then
                    THROW_HARD('Inconsistent dimensions')
                endif
                if( .not.is_equal(o%get('smpd'),o_src%get('smpd')) )then
                    THROW_HARD('Inconsistent sampling distance')
                endif
                if( o_src%isthere('movie') )then
                    call o_src%getter('movie',movfname)
                    if( movfname(1:1).ne.'/' ) movfname = trim(src_path)//'/'//trim(movfname)
                    call make_relativepath(CWD_GLOB, movfname, relmovname)
                    if( file_exists(relmovname) )then
                        call o_src%set('movie', relmovname)
                    else
                        THROW_WARN('Movie could not be substituted: '//trim(movfname))
                    endif
                endif
                if( o_src%isthere('intg') )then
                    call o_src%getter('intg',micfname)
                    if( micfname(1:1).ne.'/' ) micfname = trim(src_path)//'/'//trim(micfname)
                    call make_relativepath(CWD_GLOB, micfname, relmicname)
                    if( file_exists(relmicname) )then
                        call o_src%set('movie', relmicname)
                    else
                        call o_src%print_ori
                        write(logfhandle,*) trim(relmicname)
                        THROW_HARD('Missing micrograph!')
                    endif
                endif
                call self%os_mic%set_ori(imic,o_src)
            enddo
        case DEFAULT
            THROW_HARD('Invalid ORITYPE!')
        end select
        call o%kill
        call o_src%kill
    end subroutine replace_project

    ! projects are assumed read in
    subroutine merge_stream_projects( self, self2merge )
        class(sp_project), intent(inout) :: self, self2merge
        type(ctfparams)                  :: ctfvars, ctfvars2merge
        type(ori)                        :: o
        character(len=:),    allocatable :: output_dir, output_dir_ctf_estimate, output_dir_picker, fname
        character(len=:),    allocatable :: output_dir_motion_correct, output_dir_extract
        integer                          :: nptcls, nptcls2merge, fromp, top, top_prev, iptcl, stkind
        integer                          :: imic, nmics, nmics2merge, cnt, istk, nstks, nstks2merge, stkind_prev
        logical                          :: l_pick, l_mc, l_ctf, l_extr, err
        output_dir = PATH_HERE
        ! determines what has been done & consistency
        if( self%os_mic%get_noris() == 0 ) THROW_HARD('Source project is empty!')
        if( self2merge%os_mic%get_noris() == 0 ) THROW_HARD('Target project is empty!')
        ! motion correction
        l_mc = self%os_mic%isthere('intg') .and. self2merge%os_mic%isthere('intg')
        if( .not.l_mc )then
            THROW_HARD('One project does not contain micrographs!')
        endif
        ! ctf estimation
        ctfvars       = self%os_mic%get_ctfvars(1)
        ctfvars2merge = self2merge%os_mic%get_ctfvars(1)
        if( ctfvars%ctfflag /= ctfvars2merge%ctfflag )then
            THROW_HARD('Inconsistent CTF states between projects!')
        endif
        if( .not.is_equal(ctfvars%smpd,ctfvars2merge%smpd) )then
            THROW_HARD('Inconsistent SMPD between projects!')
        endif
        l_ctf = ctfvars%ctfflag /= CTFFLAG_NO
        ! picking
        l_pick = self%os_mic%isthere('boxfile') .or. self2merge%os_mic%isthere('boxfile')
        ! particles extraction
        l_extr = (self%os_stk%get_noris() > 0) .and. (self2merge%os_stk%get_noris() > 0)
        if( .not.l_extr )then
            THROW_WARN('Inconsistent particles extraction between projects! stacks will not be imported')
        endif
        ! folders
        output_dir_motion_correct = filepath(trim(output_dir), trim(DIR_MOTION_CORRECT))
        call simple_mkdir(output_dir_motion_correct,errmsg="commander_stream_wflows :: exec_merge_stream;  ")
        if( l_ctf )then
            output_dir_ctf_estimate = filepath(trim(output_dir), trim(DIR_CTF_ESTIMATE))
            call simple_mkdir(output_dir_ctf_estimate,errmsg="commander_stream_wflows :: exec_merge_stream;  ")
        endif
        if( l_pick )then
            output_dir_picker = filepath(trim(output_dir), trim(DIR_PICKER))
            call simple_mkdir(output_dir_picker,errmsg="commander_stream_wflows :: exec_merge_stream;  ")
        endif
        if( l_extr )then
            output_dir_extract = filepath(trim(output_dir), trim(DIR_EXTRACT))
            call simple_mkdir(output_dir_extract,errmsg="commander_stream_wflows :: exec_merge_stream;  ")
        endif
        ! micrographs segment
        nmics       = self%os_mic%get_noris()
        nmics2merge = self2merge%os_mic%get_noris()
        do imic = 1,nmics
            call self%os_mic%get_ori(imic,o)
            call movefile2folder('intg', output_dir_motion_correct, o, err)
            if( err )then
                call o%getter('intg',fname)
                THROW_HARD('Micrograph '//trim(fname)//' could not be found')
            endif
            call movefile2folder('thumb',  output_dir_motion_correct, o, err)
            call movefile2folder('forctf', output_dir_motion_correct, o, err)
            call movefile2folder('mc_starfile', output_dir_motion_correct, o, err)
            call movefile2folder('mceps',  output_dir_motion_correct, o, err)
            if( l_ctf )then
                call movefile2folder('ctfjpg', output_dir_ctf_estimate, o, err)
                call movefile2folder('ctfdoc', output_dir_ctf_estimate, o, err)
            endif
            if( l_pick )then
                call movefile2folder('boxfile',  output_dir_picker, o, err)
            endif
            call self%os_mic%set_ori(imic,o)
        enddo
        call self%os_mic%reallocate(nmics+nmics2merge)
        cnt    = nmics
        do imic = 1,nmics2merge
            cnt = cnt + 1
            call self2merge%os_mic%get_ori(imic,o)
            call movefile2folder('intg', output_dir_motion_correct, o, err)
            if( err )then
                call o%getter('intg',fname)
                THROW_HARD('Micrograph '//trim(fname)//' could not be found')
            endif
            call movefile2folder('thumb',  output_dir_motion_correct, o, err)
            call movefile2folder('forctf', output_dir_motion_correct, o, err)
            call movefile2folder('mc_starfile', output_dir_motion_correct, o, err)
            call movefile2folder('mceps',  output_dir_motion_correct, o, err)
            if( l_ctf )then
                call movefile2folder('ctfjpg', output_dir_ctf_estimate, o, err)
                call movefile2folder('ctfdoc', output_dir_ctf_estimate, o, err)
            endif
            if( l_pick )then
                call movefile2folder('boxfile',  output_dir_picker, o, err)
            endif
            call self%os_mic%set_ori(cnt,o)
        enddo
        if( l_extr )then
            ! stacks
            nstks        = self%os_stk%get_noris()
            nstks2merge  = self2merge%os_stk%get_noris()
            nptcls       = self%os_ptcl2D%get_noris()
            nptcls2merge = self2merge%os_ptcl2D%get_noris()
            do istk = 1,nstks
                call self%os_stk%get_ori(istk,o)
                call movefile2folder('stk', output_dir_extract, o, err)
                if( err )then
                    call o%getter('stk',fname)
                    THROW_HARD('Stack '//trim(fname)//' could not be found')
                endif
                call self%os_stk%set_ori(istk,o)
            enddo
            call self%os_stk%reallocate(nstks+nstks2merge)
            top_prev = nint(self%os_stk%get(nstks,'top'))
            cnt      = nstks
            do istk = 1,nstks2merge
                cnt = cnt + 1
                call self2merge%os_stk%get_ori(istk,o)
                fromp = nint(o%get('fromp')) + top_prev
                top   = nint(o%get('top'))   + top_prev
                call o%set('fromp',real(fromp))
                call o%set('top',  real(top))
                call movefile2folder('stk', output_dir_extract, o, err)
                if( err )then
                    call o%getter('stk',fname)
                    THROW_HARD('Stack '//trim(fname)//' could not be found')
                endif
                call self%os_stk%set_ori(cnt,o)
            enddo
            ! particles
            call self%os_ptcl2D%reallocate(nptcls+nptcls2merge)
            call self%os_ptcl3D%reallocate(nptcls+nptcls2merge)
            stkind_prev = nint(self%os_ptcl2D%get(nptcls,'stkind'))
            cnt    = nptcls
            do iptcl = 1,nptcls2merge
                cnt = cnt + 1
                call self2merge%os_ptcl2D%get_ori(iptcl,o)
                stkind = nint(o%get('stkind')) + stkind_prev
                call o%set('stkind',real(stkind))
                call self%os_ptcl2D%set_ori(cnt,o)
                call self2merge%os_ptcl3D%get_ori(iptcl,o)
                call o%set('stkind',real(stkind))
                call self%os_ptcl3D%set_ori(cnt,o)
            enddo
        endif
        ! the end
        call self2merge%kill
        call self%write
        call o%kill

        contains

            subroutine movefile2folder(key, folder, o, err)
                character(len=*), intent(in)    :: key, folder
                class(ori),       intent(inout) :: o
                logical,          intent(out)   :: err
                character(len=:), allocatable :: src
                character(len=LONGSTRLEN) :: dest,reldest
                integer :: iostat
                err = .false.
                if( .not.o%isthere(key) )then
                    err = .true.
                    return
                endif
                call o%getter(key,src)
                if( .not.file_exists(src) )then
                    err = .true.
                    return
                endif
                dest   = trim(folder)//'/'//basename(src)
                iostat = rename(src,dest)
                if( iostat /= 0 )then
                    THROW_WARN('Ignoring '//trim(src))
                    return
                endif
                call make_relativepath(CWD_GLOB,dest,reldest)
                call o%set(key,reldest)
            end subroutine movefile2folder

    end subroutine merge_stream_projects

    ! report state selection to os_stk & os_ptcl2D/3D
    subroutine report_state2stk( self, states )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: states(:)
        integer :: iptcl, noris_ptcl3D, noris_ptcl2D, istk, fromp, top, nstks, nptcls
        nstks = self%get_nstks()
        if( nstks == 0 )then
            THROW_WARN('empty STK field. Nothing to do; report_state2stk')
            return
        endif
        if(size(states) /= nstks)then
            THROW_WARN('Inconsistent # number of states & stacks; report_state2stk')
            return
        endif
        ! update stacks
        do istk=1,nstks
            call self%os_stk%set(istk, 'state', real(states(istk)))
        enddo
        ! ensure ptcl fields congruent
        noris_ptcl2D = self%os_ptcl2D%get_noris()
        noris_ptcl3D = self%os_ptcl3D%get_noris()
        if( noris_ptcl3D /= noris_ptcl2D )then
            THROW_HARD('Inconsistent # number of 2D/3D particles; report_state2stk')
        else
            do istk=1,nstks
                fromp  = nint(self%os_stk%get(istk,'fromp'))
                top    = nint(self%os_stk%get(istk,'top'))
                nptcls = nint(self%os_stk%get(istk,'nptcls'))
                if(top-fromp+1 /= nptcls)then
                    call self%os_stk%print_(istk)
                    THROW_HARD('Incorrect # number of particles in stack '//int2str(istk)//'; report_state2stk')
                endif
                if( states(istk) > 0 )then
                    ! preserve existing states
                    cycle
                else
                    ! de-select
                    do iptcl=fromp,top
                        call self%os_ptcl2D%set(iptcl, 'state', 0.)
                        call self%os_ptcl3D%set(iptcl, 'state', 0.)
                    enddo
                endif
            enddo
        endif
    end subroutine report_state2stk

    subroutine set_boxcoords( self, iptcl, coords )
        class(sp_project), target, intent(inout) :: self
        integer,                   intent(in)    :: iptcl, coords(2)
        call self%os_ptcl2D%set(iptcl,'xpos',real(coords(1)))
        call self%os_ptcl2D%set(iptcl,'ypos',real(coords(2)))
    end subroutine set_boxcoords

    ! printers

    subroutine print_info( self, fname )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        character(len=:),      allocatable :: projfile, record
        type(binoris_seginfo), allocatable :: hinfo(:)
        integer,               allocatable :: seginds(:)
        integer :: i
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('file format of: '//trim(fname)//' not supported; sp_project :: print_info')
        endif
        projfile = trim(fname)
        if( file_exists(projfile) )then
            call self%bos%open(projfile)
            call self%bos%get_segments_info(seginds, hinfo)
            if( allocated(hinfo) )then
                do i = 1,size(hinfo)
                    call self%bos%read_record(seginds(i), hinfo(i)%first_data_byte, record)
                    write(logfhandle,'(a)') format_str(format_str('SEGMENT '//int2str_pad(seginds(i),2)//' of type: '//segment2oritype(seginds(i)), C_BOLD), C_UNDERLINED)
                    write(logfhandle,'(a)') format_str(segment2info(seginds(i), hinfo(i)%n_records), C_ITALIC)
                    write(logfhandle,'(a)') format_str('first record:', C_BOLD)//' '//trim(record)
                end do
            endif
            call self%bos%close
        else
            THROW_HARD('projfile: '//trim(projfile)//' nonexistent; print_info')
        endif
    end subroutine print_info

    ! readers

    subroutine read( self, fname, wthreads )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        logical, optional, intent(in)    :: wthreads
        character(len=:), allocatable :: projfile
        integer :: isegment
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('format of: '//trim(fname)//' not supported; read')
        endif
        projfile = trim(fname)
        if( .not. file_exists(trim(projfile)) )then
            THROW_HARD('file: '// trim(projfile)//' does not exist; read')
        endif
        call self%bos%open(projfile)
        do isegment=1,self%bos%get_n_segments()
            call self%segreader(isegment, wthreads=wthreads)
        end do
        call self%bos%close
    end subroutine read

    subroutine read_non_data_segments( self, fname )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        character(len=:), allocatable :: projfile
        integer :: iseg
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('file format of: '//trim(fname)//' not supported; read_non_data_segments')
        endif
        projfile = trim(fname)
        if( file_exists(projfile) )then
            call self%bos%open(projfile)
            do iseg=11,MAXN_OS_SEG
                call self%segreader(iseg)
            end do
            call self%bos%close
        else
            THROW_HARD('projfile: '//trim(projfile)//' nonexistent; read_non_data_segments')
        endif
    end subroutine read_non_data_segments

    subroutine read_ctfparams_state_eo( self, fname )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        integer :: isegment
        if( .not. file_exists(trim(fname)) )then
            THROW_HARD('inputted file: '//trim(fname)//' does not exist; read_ctfparams_state_eo')
        endif
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('file format of: '//trim(fname)//' not supported; read_ctfparams_state_eo')
        endif
        call self%bos%open(fname)
        do isegment=1,self%bos%get_n_segments()
            call self%segreader(isegment, only_ctfparams_state_eo=.true.)
        end do
        call self%bos%close
    end subroutine read_ctfparams_state_eo

    subroutine read_mic_stk_ptcl2D_segments( self, fname, wthreads)
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        logical, optional, intent(in)    :: wthreads
        character(len=:), allocatable :: projfile
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('format of: '//trim(fname)//' not supported; read')
        endif
        projfile = trim(fname)
        if( .not. file_exists(trim(projfile)) )then
            THROW_HARD('file: '// trim(projfile)//' does not exist; read')
        endif
        call self%bos%open(projfile)
        call self%segreader(MIC_SEG, wthreads=wthreads)
        call self%segreader(STK_SEG, wthreads=wthreads)
        call self%segreader(PTCL2D_SEG, wthreads=wthreads)
        call self%bos%close
    end subroutine read_mic_stk_ptcl2D_segments

    subroutine read_segment( self, oritype, fname, fromto, wthreads )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: oritype
        character(len=*),  intent(in)    :: fname
        integer, optional, intent(in)    :: fromto(2)
        logical, optional, intent(in)    :: wthreads
        integer :: isegment
        if( .not. file_exists(trim(fname)) )then
            THROW_HARD('inputted file: '//trim(fname)//' does not exist; read_segment')
        endif
        select case(fname2format(fname))
            case('O')
                ! *.simple project file
                isegment = oritype2segment(oritype)
                call self%bos%open(fname)
                call self%segreader(isegment,fromto=fromto,wthreads=wthreads)
                call self%bos%close
            case('T')
                ! *.txt plain text ori file
                select case(trim(oritype))
                    case('mic')
                        call self%os_mic%read(fname)
                    case('stk')
                        call self%os_stk%read(fname)
                    case('ptcl2D')
                        call self%os_ptcl2D%read(fname, fromto=fromto)
                    case('cls2D')
                        call self%os_cls2D%read(fname)
                    case('cls3D')
                        call self%os_cls3D%read(fname,  fromto=fromto)
                    case('ptcl3D')
                        call self%os_ptcl3D%read(fname, fromto=fromto)
                    case('out')
                        call self%os_out%read(fname)
                    case('optics')
                        call self%os_optics%read(fname)
                    case('projinfo')
                        call self%projinfo%read(fname)
                    case('jobproc')
                        call self%jobproc%read(fname)
                    case('compenv')
                        call self%compenv%read(fname)
                    case DEFAULT
                        THROW_HARD('unsupported oritype flag; read_segment')
                end select
            case DEFAULT
                THROW_HARD('file format of: '//trim(fname)//' not supported; read_segment')
        end select
    end subroutine read_segment

    subroutine segreader( self, isegment, fromto, only_ctfparams_state_eo, wthreads )
        class(sp_project),          intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        integer, optional,          intent(in)    :: fromto(2)
        logical, optional,          intent(in)    :: only_ctfparams_state_eo, wthreads
        integer :: n
        n = self%bos%get_n_records(isegment)
        select case(isegment)
            case(MIC_SEG)
                call self%os_mic%new(n,    is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%os_mic,    wthreads=wthreads)
            case(STK_SEG)
                call self%os_stk%new(n,    is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%os_stk,    only_ctfparams_state_eo=only_ctfparams_state_eo, wthreads=wthreads)
            case(PTCL2D_SEG)
                call self%os_ptcl2D%new(n, is_ptcl=.true.)
                call self%bos%read_segment(isegment, self%os_ptcl2D, only_ctfparams_state_eo=only_ctfparams_state_eo, wthreads=wthreads)
            case(CLS2D_SEG)
                call self%os_cls2D%new(n,  is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%os_cls2D)
            case(CLS3D_SEG)
                call self%os_cls3D%new(n,  is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%os_cls3D)
            case(PTCL3D_SEG)
                call self%os_ptcl3D%new(n, is_ptcl=.true.)
                call self%bos%read_segment(isegment, self%os_ptcl3D, only_ctfparams_state_eo=only_ctfparams_state_eo, wthreads=wthreads)
            case(OUT_SEG)
                call self%os_out%new(n,    is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%os_out)
            case(OPTICS_SEG)
                call self%os_optics%new(n,    is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%os_optics)
            case(PROJINFO_SEG)
                call self%projinfo%new(n,  is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%projinfo)
            case(JOBPROC_SEG)
                call self%jobproc%new(n,   is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%jobproc)
            case(COMPENV_SEG)
                call self%compenv%new(n,   is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%compenv)
        end select
    end subroutine segreader


    subroutine read_segments_info( self, fname, seginds, seginfos )
        class(sp_project),                  intent(inout) :: self
        character(len=*),                   intent(in)  :: fname
        type(binoris_seginfo), allocatable, intent(out) :: seginfos(:)
        integer,               allocatable, intent(out) :: seginds(:)
        character(len=:), allocatable :: projfile
        integer :: i
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('file format of: '//trim(fname)//' not supported; sp_project :: read_segments_info')
        endif
        projfile = trim(fname)
        if( file_exists(projfile) )then
            call self%bos%open(projfile)
            call self%bos%get_segments_info(seginds, seginfos)
            call self%bos%close
        else
            THROW_HARD('projfile: '//trim(projfile)//' nonexistent; print_info')
        endif
    end subroutine read_segments_info

    !>  Convenience funtion for checking # of movies/mics, stacks and particles
    subroutine read_data_info( self, fname, nmics, nstks, nptcls )
        class(sp_project),   intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        integer,             intent(out)   :: nmics, nstks, nptcls
        type(binoris_seginfo), allocatable :: seginfos(:)
        integer,               allocatable :: seginds(:)
        integer :: i,n2D, n3D
        call self%read_segments_info(fname, seginds, seginfos)
        nmics  = 0
        nstks  = 0
        n2D    = 0
        n3D    = 0
        do i = 1,size(seginds)
            select case(seginds(i))
            case(MIC_SEG)
                nmics = int(seginfos(i)%n_records)
            case(STK_SEG)
                nstks = int(seginfos(i)%n_records)
            case(PTCL2D_SEG)
                n2D = int(seginfos(i)%n_records)
            case(PTCL3D_SEG)
                n3D = int(seginfos(i)%n_records)
            end select
        enddo
        ! if( n2D /= n3D )then
        !     THROW_WARN('Inconsistent # of particles in the 2D/3D segments; read_data_info: '//trim(fname))
        ! endif
        nptcls = n2D
    end subroutine read_data_info

    ! writers

    subroutine write( self, fname, fromto, isegment )
        class(sp_project),                    intent(inout) :: self
        character(len=*),           optional, intent(in)    :: fname
        integer,                    optional, intent(in)    :: fromto(2)
        integer(kind(ENUM_ORISEG)), optional, intent(in)    :: isegment
        character(len=:), allocatable :: projfile
        integer(kind(ENUM_ORISEG))    :: iseg
        if( present(fname) )then
            if( fname2format(fname) .ne. 'O' )then
                THROW_HARD('file format of: '//trim(fname)//' not supported; write')
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
        integer(kind(ENUM_ORISEG)) :: iseg
        if( present(fname) )then
            if( fname2format(fname) .ne. 'O' )then
                THROW_HARD('file format of: '//trim(fname)//' not supported; sp_project :: write')
            endif
            projfile = trim(fname)
        else
            call self%projinfo%getter(1, 'projfile', projfile)
        endif
        if( file_exists(projfile) )then
            iseg = oritype2segment(oritype)
            call self%bos%open(projfile)
            call self%segwriter_inside(iseg, fromto)
        else
            call self%write(fname, fromto)
        endif
        ! no need to update header (taken care of in binoris object)
        call self%bos%close
    end subroutine write_segment_inside

    subroutine write_non_data_segments( self, fname )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        character(len=:), allocatable :: projfile
        integer :: iseg
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('file format of: '//trim(fname)//' not supported; sp_project :: write_non_data_segments')
        endif
        projfile = trim(fname)
        if( file_exists(projfile) )then
            call self%bos%open(projfile)
            do iseg=11,MAXN_OS_SEG
                call self%segwriter(iseg)
            end do
            ! update header
            call self%bos%write_header
            call self%bos%close
        else
            THROW_HARD('projfile: '//trim(projfile)//' nonexistent; write_non_data_segments')
        endif
    end subroutine write_non_data_segments

    subroutine write_segment2txt( self, oritype, fname, fromto )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: oritype
        character(len=*),  intent(in)    :: fname
        integer, optional, intent(in)    :: fromto(2)
        select case(fname2format(fname))
            case('O')
                THROW_HARD('write_segment2txt is not supported for *.simple project files; write_segment2txt')
            case('T')
                ! *.txt plain text ori file
                select case(trim(oritype))
                    case('mic')
                        if( self%os_mic%get_noris() > 0 )then
                            call self%os_mic%write(fname)
                        else
                            THROW_WARN('no mic-type oris available to write; write_segment2txt')
                        endif
                    case('stk')
                        if( self%os_stk%get_noris() > 0 )then
                            call self%os_stk%write(fname)
                        else
                            THROW_WARN('no stk-type oris available to write; write_segment2txt')
                        endif
                    case('ptcl2D')
                        if( self%os_ptcl2D%get_noris() > 0 )then
                            call self%os_ptcl2D%write(fname, fromto)
                        else
                            THROW_WARN('no ptcl2D-type oris available to write; write_segment2txt')
                        endif
                    case('cls2D')
                        if( self%os_cls2D%get_noris() > 0 )then
                            call self%os_cls2D%write(fname)
                        else
                            THROW_WARN('no cls2D-type oris available to write; write_segment2txt')
                        endif
                    case('cls3D')
                        if( self%os_cls3D%get_noris() > 0 )then
                            call self%os_cls3D%write(fname,  fromto)
                        else
                            THROW_WARN('no cls3D-type oris available to write; write_segment2txt')
                        endif
                    case('ptcl3D')
                        if( self%os_ptcl3D%get_noris() > 0 )then
                            call self%os_ptcl3D%write(fname, fromto)
                        else
                            THROW_WARN('no ptcl3D-type oris available to write; write_segment2txt')
                        endif
                    case('out')
                        if( self%os_out%get_noris() > 0 )then
                            call self%os_out%write(fname)
                        else
                            THROW_WARN('no out-type oris available to write; write_segment2txt')
                        endif
                    case('optics')
                        if( self%os_optics%get_noris() > 0 )then
                            call self%os_optics%write(fname)
                        else
                            THROW_WARN('no optics-type oris available to write; write_segment2txt')
                        endif
                    case('projinfo')
                        if( self%projinfo%get_noris() > 0 )then
                            call self%projinfo%write(fname, fromto)
                        else
                            THROW_WARN('no projinfo-type oris available to write; write_segment2txt')
                        endif
                    case('jobproc')
                        if( self%jobproc%get_noris() > 0 )then
                            call self%jobproc%write(fname)
                        else
                            THROW_WARN('no jobproc-type oris available to write; write_segment2txt')
                        endif
                    case('compenv')
                        if( self%compenv%get_noris() > 0 )then
                            call self%compenv%write(fname)
                        else
                            THROW_WARN('no compenv-type oris available to write; write_segment2txt')
                        endif
                    case DEFAULT
                        THROW_HARD('unsupported oritype flag; write_segment2txt')
                end select
            case DEFAULT
                THROW_HARD('file format of: '//trim(fname)//'not supported; write_segment2txt')
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
                        write(logfhandle,'(a)') self%os_mic%ori2str(iori)
                    end do
                else
                    write(logfhandle,*) 'No mic-type oris available to print; sp_project :: print_segment'
                endif
            case('stk')
                noris = self%os_stk%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(logfhandle,'(a)') self%os_stk%ori2str(iori)
                    end do
                else
                    write(logfhandle,*) 'No stk-type oris available to print; sp_project :: print_segment'
                endif
            case('ptcl2D')
                noris = self%os_ptcl2D%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(logfhandle,'(a)') self%os_ptcl2D%ori2str(iori)
                    end do
                else
                    write(logfhandle,*) 'No ptcl2D-type oris available to print; sp_project :: print_segment'
                endif
            case('cls2D')
                noris = self%os_cls2D%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(logfhandle,'(a)') self%os_cls2D%ori2str(iori)
                    end do
                else
                    write(logfhandle,*) 'No cls2D-type oris available to print; sp_project :: print_segment'
                endif
            case('cls3D')
                noris = self%os_cls3D%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(logfhandle,'(a)') self%os_cls3D%ori2str(iori)
                    end do
                else
                    write(logfhandle,*) 'No cls3D-type oris available to print; sp_project :: print_segment'
                endif
            case('ptcl3D')
                noris = self%os_ptcl3D%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(logfhandle,'(a)') self%os_ptcl3D%ori2str(iori)
                    end do
                else
                    write(logfhandle,*) 'No ptcl3D-type oris available to print; sp_project :: print_segment'
                endif
            case('out')
                noris = self%os_out%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(logfhandle,'(a)') self%os_out%ori2str(iori)
                    end do
                else
                    write(logfhandle,*) 'No out-type oris available to print; sp_project :: print_segment'
                endif
            case('optics')
                noris = self%os_optics%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(logfhandle,'(a)') self%os_optics%ori2str(iori)
                    end do
                else
                    write(logfhandle,*) 'No optics-type oris available to print; sp_project :: print_segment'
                endif
            case('projinfo')
                noris = self%projinfo%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(logfhandle,'(a)') self%projinfo%ori2str(iori)
                    end do
                else
                    write(logfhandle,*) 'No projinfo-type oris available to print; sp_project :: print_segment'
                endif
            case('jobproc')
                noris = self%jobproc%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(logfhandle,'(a)') self%jobproc%ori2str(iori)
                    end do
                else
                    write(logfhandle,*) 'No jobproc-type oris available to print; sp_project :: print_segment'
                endif
            case('compenv')
                noris = self%compenv%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        write(logfhandle,'(a)') self%compenv%ori2str(iori)
                    end do
                else
                    write(logfhandle,*) 'No compenv-type oris available to print; sp_project :: print_segment'
                endif
            case DEFAULT
                THROW_HARD('unsupported oritype flag; print_segment')
        end select
    end subroutine print_segment

    subroutine segwriter( self, isegment, fromto )
        class(sp_project), intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
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
            case(OPTICS_SEG)
                call self%bos%write_segment(isegment, self%os_optics)
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
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
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
            case(OPTICS_SEG)
                call self%bos%write_segment_inside(isegment, self%os_optics)
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
        call self%os_optics%kill
        call self%projinfo%kill
        call self%jobproc%kill
        call self%compenv%kill
    end subroutine kill

    ! private supporting subroutines / functions

    integer(kind(ENUM_ORISEG)) function oritype2segment( oritype )
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
            case('optics')
                oritype2segment = OPTICS_SEG
            case('projinfo')
                oritype2segment = PROJINFO_SEG
            case('jobproc')
                oritype2segment = JOBPROC_SEG
            case('compenv')
                oritype2segment = COMPENV_SEG
            case DEFAULT
                THROW_HARD('unsupported oritype flag; oritype2segment')
        end select
    end function oritype2segment

    function segment2oritype( iseg ) result( oritype )
        integer(kind(ENUM_ORISEG)), intent(in) :: iseg
        character(len=:), allocatable :: oritype
        select case(iseg)
            case(MIC_SEG)
                oritype = 'mic'
            case(STK_SEG)
                oritype = 'stk'
            case(PTCL2D_SEG)
                oritype = 'ptcl2D'
            case(CLS2D_SEG)
                oritype = 'cls2D'
            case(CLS3D_SEG)
                oritype = 'cls3D'
            case(PTCL3D_SEG)
                oritype = 'ptcl3D'
            case(OUT_SEG)
                oritype = 'out'
            case(OPTICS_SEG)
                oritype = 'optics'
            case(PROJINFO_SEG)
                oritype = 'projinfo'
            case(JOBPROC_SEG)
                oritype = 'jobproc'
            case(COMPENV_SEG)
                oritype = 'compenv'
            case DEFAULT
                write(logfhandle,*) 'iseg: ', iseg
                THROW_HARD('unsupported segment of kind(ENUM_ORISEG); segment2oritype')
        end select
    end function segment2oritype

    function segment2info( iseg, n_records ) result( info )
        integer(kind(ENUM_ORISEG)), intent(in) :: iseg
        integer(kind=8),            intent(in) :: n_records
        character(len=:), allocatable :: info, nrecs_str
        nrecs_str = int2str(int(n_records,kind=4))
        info = nrecs_str//' record(s) of '
        select case(iseg)
            case(MIC_SEG)
                info = info//'movie and micrograph (integrated movie) info, one per movie/micrograph'
            case(STK_SEG)
                info = info//'stack (extracted particles) info, one per stack of particles'
            case(PTCL2D_SEG)
                info = info//'2D information generated by cluster2D, one per particle'
            case(CLS2D_SEG)
                info = info//'data generated by cluster2D, one per 2D cluster'
            case(CLS3D_SEG)
                info = info//'3D information for class averages, one per class'
            case(PTCL3D_SEG)
                info = info//'3D information, one per particle'
            case(OUT_SEG)
                info = info//'crtitical project outputs: class averages, 3D volumes, FSC/FRC files etc.'
            case(OPTICS_SEG)
                info = info//'optics group information.'
            case(PROJINFO_SEG)
                info = info//'information about the project, project name etc.'
            case(JOBPROC_SEG)
                info = info//'all command-lines executed throughout the project'
            case(COMPENV_SEG)
                info = info//'computing environment specifications, queue system, memory per job etc.'
            case DEFAULT
                write(logfhandle,*) 'iseg: ', iseg
                THROW_HARD('unsupported segment of kind(ENUM_ORISEG); segment2oritype')
        end select
    end function segment2info

end module simple_sp_project
