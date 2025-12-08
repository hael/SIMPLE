submodule(simple_sp_project) simple_sp_project_core
! include 'simple_lib.f08'
implicit none
#include "simple_local_flags.inc"
contains

    module subroutine new_seg_with_ptr( self, n, oritype, os_ptr )
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

    module subroutine kill( self )
        class(sp_project), intent(inout) :: self
        call self%os_mic%kill
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
        call self%bos%close
    end subroutine kill

    module subroutine assign( self_out, self_in )
        class(sp_project), intent(inout) :: self_out
        class(sp_project), intent(in)    :: self_in
        call self_out%copy(self_in)
    end subroutine assign

    module subroutine copy( self_out, self_in )
        class(sp_project), intent(inout) :: self_out
        class(sp_project), intent(in)    :: self_in
        call self_out%os_mic%copy(self_in%os_mic)
        call self_out%os_stk%copy(self_in%os_stk)
        call self_out%os_ptcl2D%copy(self_in%os_ptcl2D)
        call self_out%os_cls2D%copy(self_in%os_cls2D)
        call self_out%os_ptcl3D%copy(self_in%os_ptcl3D)
        call self_out%os_cls3D%copy(self_in%os_cls3D)
        call self_out%os_out%copy(self_in%os_out)
        call self_out%os_optics%copy(self_in%os_optics)
        call self_out%projinfo%copy(self_in%projinfo)
        call self_out%jobproc%copy(self_in%jobproc)
        call self_out%compenv%copy(self_in%compenv)
    end subroutine copy

    ! field updaters

    module subroutine update_projinfo_1( self, cline )
        class(sp_project), intent(inout) :: self
        class(cmdline),    intent(in)    :: cline
        type(string) :: projname, projfile, cwd
        if( self%projinfo%get_noris() == 1 )then
            ! no need to construct field
        else
            call self%projinfo%new(1, is_ptcl=.false.)
        endif
        ! projname & profile
        if( self%projinfo%isthere('projname').and.cline%defined('projname') )then
            projname = cline%get_carg('projname')
            call self%projinfo%set(1, 'projname', projname%to_char())
            call self%projinfo%set(1, 'projfile', projname%to_char()//'.simple')
        else
            if( .not. cline%defined('projname') .and. .not. cline%defined('projfile') )then
                THROW_HARD('the project needs a name, inputted via projname or projfile; update_projinfo')
            endif
            if( cline%defined('projfile') )then
                projfile = cline%get_carg('projfile')
                select case(fname2format(projfile))
                    case('O')
                        call self%projinfo%set(1, 'projfile', projfile%to_char())
                    case DEFAULT
                        THROW_HARD('unsupported format of projfile: '//projfile%to_char()//'; update_projinfo')
                end select
                projname = get_fbody(projfile, string('simple'))
                call self%projinfo%set(1, 'projname', projname%to_char())
            endif
            if( cline%defined('projname') )then
                projname = cline%get_carg('projname')
                call self%projinfo%set(1, 'projname', projname%to_char())
                call self%projinfo%set(1, 'projfile', projname%to_char()//'.simple')
            endif
        endif
        ! it is assumed that the project is created in the root "project directory", i.e. stash cwd
        call simple_getcwd(cwd)
        call self%projinfo%set(1, 'cwd', cwd%to_char())
    end subroutine update_projinfo_1

    module subroutine update_projinfo_2( self, projfile )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: projfile
        type(string) :: projname, cwd
        if( self%projinfo%get_noris() == 1 )then
            ! no need to construct field
        else
            call self%projinfo%new(1, is_ptcl=.false.)
        endif
        ! projfile & projname
        select case(fname2format(projfile))
            case('O')
                call self%projinfo%set(1, 'projfile', projfile%to_char() )
            case DEFAULT
                THROW_HARD('unsupported format of projfile: '//projfile%to_char()//'; update_projinfo_2')
        end select
        projname = get_fbody(projfile, string('simple'))
        call self%projinfo%set(1, 'projname', projname%to_char())
        ! it is assumed that the project is created in the root "project directory", i.e. stash cwd
        call simple_getcwd(cwd)
        call self%projinfo%set(1, 'cwd', cwd%to_char())
    end subroutine update_projinfo_2

    module subroutine update_compenv( self, cline )
        class(sp_project), intent(inout) :: self
        class(cmdline),    intent(in)    :: cline
        type(string) :: env_var
        type(string) :: projname, qsnam
        integer      :: iostat
        if( self%compenv%get_noris() == 1 )then
            ! no need to construct field
        else
            call self%compenv%new(1, is_ptcl=.false.)
        endif
        ! compenv has to be filled as strings as it is used as a string only dictionary
        ! get from environment
        env_var = simple_getenv('SIMPLE_PATH', iostat)
        if( iostat /= 0 )then
            write(logfhandle,*) 'ERROR! SIMPLE_PATH is not defined in your shell environment!'
            write(logfhandle,*) 'Please refer to installation documentation for correct system configuration'
            stop
        else
            call self%compenv%set(1, 'simple_path', env_var)
        endif
        if( cline%defined('qsys_name') )then
            qsnam = cline%get_carg('qsys_name')
            call self%compenv%set(1, 'qsys_name', qsnam%to_char())
            iostat = 0
        else
            env_var = simple_getenv('SIMPLE_QSYS', iostat)
            if( iostat == 0 ) call self%compenv%set(1, 'qsys_name', env_var)
        endif
        if( iostat /= 0 ) THROW_HARD('SIMPLE_QSYS is not defined in your environment; update_compenv')
        env_var = simple_getenv('SIMPLE_EMAIL', iostat)
        if( iostat/=0 ) env_var = 'my.name@uni.edu'
        ! get from command line
        call self%compenv%set(1, 'user_email', env_var)
        if( cline%defined('user_email') )then
            call self%compenv%set(1, 'user_email', cline%get_carg('user_email'))
        else
            call self%compenv%set(1, 'user_email', env_var)
        endif
        if( cline%defined('time_per_image') )then
            call self%compenv%set(1, 'time_per_image', cline%get_iarg('time_per_image'))
        else
            if( .not. self%compenv%isthere('time_per_image') )then
                call self%compenv%set(1, 'time_per_image', TIME_PER_IMAGE_DEFAULT)
            endif
        endif
        if( cline%defined('walltime') )then
            call self%compenv%set(1, 'walltime', cline%get_iarg('walltime'))
        else
            if( .not. self%compenv%isthere('walltime') )then
                call self%compenv%set(1, 'walltime', WALLTIME_DEFAULT)
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
        else
            env_var = simple_getenv('SIMPLE_QSYS_PARTITION', iostat)
            if( iostat == 0 ) call self%compenv%set(1, 'qsys_partition', env_var)
        endif
        if( cline%defined('qsys_qos') )then
            call self%compenv%set(1, 'qsys_qos', cline%get_carg('qsys_qos'))
        endif
        if( cline%defined('qsys_reservation') )then
            call self%compenv%set(1, 'qsys_reservation', cline%get_carg('qsys_reservation'))
        endif
        if( .not. self%compenv%isthere('job_name') )then
            call self%projinfo%getter(1, 'projname', projname)
            call self%compenv%set(1, 'job_name', 'simple_'//projname%to_char() )
        endif
        if( cline%defined('job_memory_per_task') )then
            call self%compenv%set(1, 'job_memory_per_task', cline%get_iarg('job_memory_per_task') )
        else
            if( .not. self%compenv%isthere('job_memory_per_task') )then
                call self%compenv%set(1, 'job_memory_per_task', JOB_MEMORY_PER_TASK_DEFAULT)
            endif
        endif
    end subroutine update_compenv

    !>  append project2 to project1
    !   projinfo,jobproc,compenv fields from project2 are ignored
    !   cls2D/cls3D/out fields are wiped
    module subroutine append_project( self1, self2 )
        class(sp_project), intent(inout) :: self1
        class(sp_project), intent(in)    :: self2
        integer :: nmics1, nstks1, nptcls1, nmics2, nstks2, nptcls2, nogs1, nogs2
        integer :: i, iptcl, imic, istk, fromp, top, og_offset, ogid
        logical :: l_has_mics, l_has_stks, l_has_ptcls, l_has_optics
        nmics1  = self1%os_mic%get_noris()
        nmics2  = self2%os_mic%get_noris()
        nstks1  = self1%get_nstks()
        nstks2  = self2%get_nstks()
        nptcls1 = self1%get_nptcls()
        nptcls2 = self2%get_nptcls()
        nogs1   = self1%os_optics%get_noris()
        nogs2   = self2%os_optics%get_noris()
        ! sanity checks
        if( (nmics1==0) .neqv.(nmics2==0)  ) THROW_HARD('Only one mic field is populated!')
        if( (nstks1==0) .neqv.(nstks2==0)  ) THROW_HARD('Only one stk field is populated!')
        if( (nptcls1==0).neqv.(nptcls2==0) ) THROW_HARD('Only one ptcl field is populated!')
        if( (nogs1==0)  .neqv.(nogs2==0)   ) THROW_HARD('Only one optics field is populated!')
        l_has_mics   = (nmics1 > 0)  .and. (nmics2 > 0)
        l_has_stks   = (nstks1 > 0)  .and. (nstks2 > 0)
        l_has_ptcls  = (nptcls1 > 0) .and. (nptcls2 > 0)
        l_has_optics = (nogs1 > 0)   .and. (nogs2 > 0)
        if( l_has_stks .neqv. l_has_ptcls ) THROW_HARD('Missing stk/ptcl field!')
        ! micrograph field
        if( l_has_mics )then
            call self1%os_mic%append(self2%os_mic)
            write(logfhandle,'(A,I8)')'>>> CURRENT # OF MOVIES/MICROGRAPHS: ',nmics1+nmics2
        endif
        ! stacks & particles
        if( l_has_stks )then
            ! stack
            call self1%os_stk%append(self2%os_stk)
            write(logfhandle,'(A,I8)')'>>> CURRENT # OF STACKS:             ',nstks1+nstks2
            ! particles
            call self1%os_ptcl2D%append(self2%os_ptcl2D)
            call self1%os_ptcl3D%append(self2%os_ptcl3D)
            ! update particles/stk indices
            do istk = nstks1+1,nstks1+nstks2
                fromp = self1%os_stk%get_fromp(istk)
                top   = self1%os_stk%get_top(istk)
                fromp = fromp + nptcls1
                top   = top   + nptcls1
                call self1%os_stk%set(istk, 'fromp', fromp)
                call self1%os_stk%set(istk, 'top',   top)
                do iptcl = fromp,top
                    call self1%os_ptcl2D%set_stkind(iptcl, istk)
                    call self1%os_ptcl3D%set_stkind(iptcl, istk)
                enddo
            enddo
            write(logfhandle,'(A,I8)')'>>> CURRENT # OF PARTICLES:          ',nptcls1+nptcls2
        endif
        ! cls2D/3D
        call self1%os_cls2D%kill
        call self1%os_cls3D%kill
        ! out is wiped
        call self1%os_out%kill
        ! optic groups
        if( l_has_optics )then
            call self1%os_optics%append(self2%os_optics)
            ! determining numbering offset
            og_offset = -1
            do i = 1,nogs1
                og_offset = max(og_offset, self1%os_optics%get_int(i,'ogid'))
            end do
            if( og_offset < 1 ) THROW_HARD('Invalid optics field!')
            og_offset = max(nogs1, og_offset)
            ! updating os_optics
            do i = nogs1+1,nogs1+nogs2
                ogid = self1%os_optics%get_int(i,'ogid') + og_offset
                call self1%os_optics%set(i,'ogid',  ogid)
                call self1%os_optics%set(i,'ogname','opticsgroup'//int2str(ogid))
            enddo
            if( l_has_mics )then
                ! updating os_mic
                do imic = nmics1+1,nmics1+nmics2
                    ogid = self1%os_mic%get_int(imic,'ogid') + og_offset
                    call self1%os_mic%set(imic,'ogid', ogid)
                enddo
            endif
            if( l_has_ptcls )then
                ! updating particles
                do iptcl = nptcls1+1,nptcls1+nptcls2
                    ogid = self1%os_ptcl2D%get_int(iptcl,'ogid') + og_offset
                    call self1%os_ptcl2D%set(iptcl, 'ogid', ogid)
                    call self1%os_ptcl3D%set(iptcl, 'ogid', ogid)
                enddo
            endif
            write(logfhandle,'(A,I8)')'>>> CURRENT # OF OPTICS GROUPS:      ',nogs1+nogs2
        endif
    end subroutine append_project

    module subroutine append_job_descr2jobproc( self, exec_dir, job_descr, did_update )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: exec_dir
        class(chash),      intent(inout) :: job_descr
        logical,           intent(out)   :: did_update
        type(string)  :: edir
        type(ori)     :: o
        integer       :: njobs, ijob, ind
        character(8)  :: date
        character(10) :: time
        did_update = .true.
        njobs = self%jobproc%get_noris()
        if( njobs > 0 )then
            if( exec_dir%to_char() .ne. './')then
                do ijob=1,njobs
                    if( self%jobproc%isthere(ijob, 'exec_dir') )then
                        call self%jobproc%getter(ijob, 'exec_dir', edir)
                        if( exec_dir%has_substr(edir) )then
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
        call o%set('exec_dir', exec_dir%to_char())
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

    module subroutine replace_project( self, projfile_src, oritype )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: projfile_src
        character(len=*),  intent(in)    :: oritype
        type(string)    :: absstkname, absboxname, absmicname, absmovname, boxfname, stkfname, movfname, micfname, src_path
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
            if( self_src%os_stk%get_top(nstks) /= self%os_ptcl2D%get_noris() )then
                THROW_HARD('Inconsistent # of particles between stk & particles fields')
            endif
            if( nstks == 0 ) return
            do istk = 1,nstks
                err   = .false.
                call self_src%os_stk%get_ori(istk, o_src)
                call self%os_stk%get_ori(istk, o)
                if( o%get_int('box') /= o_src%get_int('box') )then
                    THROW_HARD('Inconsistent box size')
                endif
                if( .not.is_equal(o%get('smpd'),o_src%get('smpd')) )then
                    THROW_HARD('Inconsistent sampling distance')
                endif
                if( o_src%isthere('stk') )then
                    call o_src%getter('stk',stkfname)
                    if( stkfname%to_char([1,1]).ne.'/' ) stkfname = src_path%to_char()//'/'//stkfname%to_char()
                    absstkname = simple_abspath(stkfname, check_exists=.false.)
                    if( file_exists(absstkname) )then
                        call o_src%set('stk', absstkname)
                    else
                        err = .true.
                    endif
                else
                    err = .true.
                endif
                if( o_src%isthere('boxfile') )then
                    call o_src%getter('boxfile',boxfname)
                    if( boxfname%to_char([1,1]).ne.'/' ) boxfname = src_path%to_char()//'/'//boxfname%to_char()
                    absboxname = simple_abspath(boxfname, check_exists=.false.)
                    if( file_exists(absboxname) )then
                        call o_src%set('boxfile', absboxname)
                    else
                        err = .true.
                    endif
                endif
                if( err )then
                    call o_src%print_ori
                    write(logfhandle,*) absstkname%to_char()
                    write(logfhandle,*) absboxname%to_char()
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
                if(  o%get_int('xdim') /= o_src%get_int('xdim') .or.&
                    &o%get_int('ydim') /= o_src%get_int('ydim'))then
                    THROW_HARD('Inconsistent dimensions')
                endif
                if( .not.is_equal(o%get('smpd'),o_src%get('smpd')) )then
                    THROW_HARD('Inconsistent sampling distance')
                endif
                if( o_src%isthere('movie') )then
                    call o_src%getter('movie',movfname)
                    if( movfname%to_char([1,1]).ne.'/' ) movfname = src_path%to_char()//'/'//movfname%to_char()
                    absmovname = simple_abspath(movfname, check_exists=.false.)
                    if( file_exists(absmovname) )then
                        call o_src%set('movie', absmovname)
                    else
                        THROW_WARN('Movie could not be substituted: '//absmovname%to_char())
                    endif
                endif
                if( o_src%isthere('intg') )then
                    call o_src%getter('intg',micfname)
                    if( micfname%to_char([1,1]).ne.'/' ) micfname = src_path%to_char()//'/'//micfname%to_char()
                    absmicname = simple_abspath(micfname, check_exists=.false.)
                    if( file_exists(absmicname) )then
                        call o_src%set('movie', absmicname)
                    else
                        call o_src%print_ori
                        write(logfhandle,*) absmicname%to_char()
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

    !> for merging alignment documents from SIMPLE runs in distributed mode
    module subroutine merge_algndocs( self, nptcls, ndocs, oritype, fbody, numlen_in )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: nptcls, ndocs
        character(len=*),  intent(in)    :: oritype, fbody
        integer, optional, intent(in)    :: numlen_in
        integer,      allocatable :: parts(:,:)
        type(string)              :: fname, projfile
        type(string), allocatable :: os_strings(:)
        class(oris),      pointer :: os => null()
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
                call bos_doc%open(fname)
                n_records = bos_doc%get_n_records(isegment)
                partsz    = parts(i,2) - parts(i,1) + 1
                if( n_records /= partsz )then
                    write(logfhandle,*) 'ERROR, # records does not match expectation'
                    write(logfhandle,*) 'EXTRACTED FROM file: ', fname%to_char()
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
                call bos_doc%open(fname)
                n_records = bos_doc%get_n_records(isegment)
                partsz    = parts(i,2) - parts(i,1) + 1
                if( n_records /= partsz )then
                    write(logfhandle,*) 'ERROR, # records does not match expectation'
                    write(logfhandle,*) 'EXTRACTED FROM file: ', fname%to_char()
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
                call os%str2ori(i, os_strings(i)%to_char())
            end do
            call os_strings%kill
        endif
        if( allocated(parts) ) deallocate(parts)
        nullify(os)
        ! no need to update header (taken care of in binoris object)
        call self%bos%close
    end subroutine merge_algndocs

    !> convert a polymorphic list of project records into a sp_project instance
    !  previous mic/stk/ptcl2D,ptcl3D are wipped, other fields untouched
    module subroutine projrecords2proj( spproj, project_list )
        class(sp_project), intent(inout) :: spproj
        class(rec_list ),  intent(inout) :: project_list
        type(project_rec)  :: prec
        type(rec_iterator) :: it
        type(sp_project)   :: tmpproj
        type(string)       :: stack_name, projname, prev_projname
        integer :: iptcl, fromp, ifromp, itop, jptcl, nptcls_tot
        integer :: nrecs, nmics, nptcls, imic, micind
        logical :: has_ptcl
        call spproj%os_mic%kill
        call spproj%os_stk%kill
        call spproj%os_ptcl2D%kill
        call spproj%os_ptcl3D%kill
        nrecs      = project_list%size()
        if( nrecs == 0 ) return 
        nmics      = nrecs
        nptcls_tot = project_list%get_nptcls_tot()
        has_ptcl   = nptcls_tot > 0
        call spproj%os_mic%new(nmics,is_ptcl=.false.)
        call spproj%os_stk%new(nmics,is_ptcl=.false.)
        if( has_ptcl ) call spproj%os_ptcl2D%new(nptcls_tot,is_ptcl=.true.)
        prev_projname = ''
        jptcl = 0
        fromp = 1
        it    = project_list%begin()
        do imic = 1,nmics
            ! read individual project (up to STREAM_NMOVS_SET entries)
            ! retrieve one record from the list with the iterator
            call it%get(prec)
            projname = prec%projname
            if( projname /= prev_projname )then
                call tmpproj%kill
                call tmpproj%read_mic_stk_ptcl2D_segments(projname)
                prev_projname = projname
            endif
            ! mic
            micind = prec%micind
            call spproj%os_mic%transfer_ori(imic, tmpproj%os_mic, micind)
            ! stack
            nptcls = prec%nptcls
            if( nptcls == 0 )then
                ! move the iterator
                call it%next()
                cycle
            endif
            call spproj%os_stk%transfer_ori(imic, tmpproj%os_stk, micind)
            ! update stack path to absolute
            stack_name = spproj%get_stkname(imic)
            if( stack_name%to_char([1,1]) == '/' )then
                ! already absolute path, should always be the case
            else if( stack_name%to_char([1,3]) == '../' )then
                stack_name = simple_abspath(stack_name)
                call spproj%os_stk%set(imic, 'stk', stack_name)
            else
                THROW_HARD('Unexpected file path format for: '//stack_name%to_char())
            endif
            ! particles
            ifromp = spproj%os_stk%get_fromp(imic)
            itop   = spproj%os_stk%get_top(imic)
            do iptcl = ifromp,itop
                jptcl = jptcl+1 ! global index
                call spproj%os_ptcl2D%transfer_ori(jptcl, tmpproj%os_ptcl2D, iptcl)
                call spproj%os_ptcl2D%set_stkind(jptcl, imic)
            enddo
            call spproj%os_stk%set(imic, 'fromp', fromp)
            call spproj%os_stk%set(imic, 'top',   fromp+nptcls-1)
            fromp = fromp + nptcls
            ! move the iterator
            call it%next()
        enddo
        call tmpproj%kill
        if( has_ptcl ) spproj%os_ptcl3D = spproj%os_ptcl2D
    end subroutine projrecords2proj

    ! Getters/Setters

    module subroutine ptr2oritype( self, oritype, os_ptr )
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

    module real function get_smpd( self )
        class(sp_project), target, intent(inout) :: self
        integer :: n_os_stk, n_os_mic
        get_smpd = 0.
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk == 0 )then
            n_os_mic = self%os_mic%get_noris()
            if( n_os_mic == 0 )then
                THROW_HARD('empty os_mic field! get_smpd, aborting')
            else
                get_smpd = self%os_mic%get(1,'smpd')
            endif
        else
            get_smpd = self%os_stk%get(1,'smpd')
        endif
    end function get_smpd

    ! static for OpenMP safety
    module function get_ctfflag( self, oritype, iptcl ) result( ctfflag )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,         optional, intent(in)    :: iptcl
        character(len=STDLEN) :: ctfflag
        class(oris), pointer  :: ptcl_field
        type(string)          :: ctfflag_str
        integer :: stkind, ind_in_stk, ind
        ctfflag = 'no'
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
        ctfflag = ctfflag_str%to_char()
    end function get_ctfflag

    module integer(kind(ENUM_CTFFLAG)) function get_ctfflag_type( self, oritype, iptcl )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,         optional, intent(in)    :: iptcl
        character(len=STDLEN) :: ctfflag
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

    module integer function get_n_insegment( self, oritype )
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

    module integer function get_n_insegment_state( self, oritype, state )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,                   intent(in)    :: state
        class(oris), pointer :: pos => NULL()
        integer              :: iori
        get_n_insegment_state = 0
        select case(oritype)
            case('ptcl2D','ptcl3D')
                ! # ptcl2D = # ptcl3D
                if( self%os_ptcl2D%get_noris() /= self%os_ptcl3D%get_noris() )then
                   THROW_HARD('Inconsistent number of particles in PTCL2D/PTCL3D segments; get_n_insegment_state')
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

    ! static for OpenMP safety
    module function get_ctfparams( self, oritype, iptcl ) result( ctfvars )
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
            call self%os_stk%get_static(stkind, 'ctf', ctfflag)
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
        select case(trim(oritype))
            case('stk')
                ! nothing to do
            case DEFAULT
                ! defocus in x
                if( ptcl_field%isthere(iptcl, 'dfx') )then
                    ctfvars%dfx = ptcl_field%get_dfx(iptcl)
                else
                    call ptcl_field%print(iptcl)
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
                        call ptcl_field%print(iptcl)
                        THROW_HARD('astigmatic CTF model requires angast (angle of astigmatism) lacking in os_stk field; get_ctfparams')
                    else
                        ctfvars%angast = 0.
                    endif
                endif
        end select
        ! has phaseplate
        if( self%os_stk%isthere(stkind, 'phaseplate') )then
            call self%os_stk%get_static(stkind, 'phaseplate', phaseplate)
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

    module subroutine get_sp_oris( self, which_imgkind, os )
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

    module subroutine set_sp_oris( self, which_imgkind, os )
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

end submodule simple_sp_project_core
