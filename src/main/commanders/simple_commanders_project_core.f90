!@descr: commanders for operating on projects (spproject) and associated files, the core stuff
module simple_commanders_project_core
use simple_commanders_api
use simple_stream_communicator, only: stream_http_communicator
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_new_project
  contains
    procedure :: execute      => exec_new_project
end type commander_new_project

type, extends(commander_base) :: commander_print_project_info
  contains
    procedure :: execute      => exec_print_project_info
end type commander_print_project_info

type, extends(commander_base) :: commander_print_project_vals
  contains
    procedure :: execute      => exec_print_project_vals
end type commander_print_project_vals

type, extends(commander_base) :: commander_print_project_field
  contains
    procedure :: execute      => exec_print_project_field
end type commander_print_project_field

type, extends(commander_base) :: commander_update_project
  contains
    procedure :: execute      => exec_update_project
end type commander_update_project

type, extends(commander_base) :: commander_merge_projects
  contains
    procedure :: execute      => exec_merge_projects
end type commander_merge_projects

type, extends(commander_base) :: commander_replace_project_field
  contains
    procedure :: execute      => exec_replace_project_field
end type commander_replace_project_field

type, extends(commander_base) :: commander_concatenate_projects
  contains
    procedure :: execute      => exec_concatenate_projects
end type commander_concatenate_projects

type, extends(commander_base) :: commander_aggregate_chunks
  contains
    procedure :: execute      => exec_aggregate_chunks
end type commander_aggregate_chunks

type, extends(commander_base) :: commander_extract_subproj
contains
    procedure :: execute      => exec_extract_subproj
end type commander_extract_subproj

type, extends(commander_base) :: commander_selection
  contains
    procedure :: execute      => exec_selection
end type commander_selection

contains

    subroutine exec_new_project( self, cline )
        class(commander_new_project), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        call cline%set('mkdir', 'no')
        call params%new(cline)
        if( cline%defined('projfile') )then
            call spproj%read(params%projfile)
        endif
        if( cline%defined('projname') .and. cline%defined('dir') )then
            if( params%projname%substr_ind('.') > 0 ) THROW_HARD('No punctuation allowed in projname')
            if( .not. file_exists(params%dir) )then
                write(logfhandle,*) 'input project directory (dir): ', params%dir%to_char(), ' does not exist'
                THROW_HARD('ABORTING... exec_new_project')
            endif
            if( file_exists(params%dir // '/' // params%projname%to_char() // ".simple") )then
                write(logfhandle,*) 'input project file : ', params%dir%to_char() // '/' // params%projname%to_char() //".simple", ' already exists'
                THROW_HARD('ABORTING... exec_new_project')
            endif
            ! change to project directory
            call simple_chdir(params%dir)
            call cline%set('projname', params%projname)
        else if( cline%defined('projname') )then
            if( params%projname%substr_ind('.') > 0 ) THROW_HARD('No punctuation allowed in projname')
            if( file_exists(string(PATH_HERE)//params%projname) )then
                write(logfhandle,*) 'project directory: ', params%projname%to_char(), ' already exists in cwd: ', params%cwd%to_char()
                write(logfhandle,*) 'If you intent to overwrite the existing file, please remove it and re-run new_project'
                THROW_HARD('ABORTING... exec_new_project')
            endif
            ! make project directory
            call simple_mkdir(params%projname)
            ! change to project directory
            call simple_chdir(params%projname)
        else if( cline%defined('dir') )then
            params%dir = simple_abspath(params%dir)
            if( .not. file_exists(params%dir) )then
                write(logfhandle,*) 'input project directory (dir): ', params%dir%to_char(), ' does not exist'
                THROW_HARD('ABORTING... exec_new_project')
            endif
            call cline%set('projname', basename(params%dir))
            ! change to project directory
            call simple_chdir(params%dir)
        else
            THROW_HARD('neither projname nor dir defined on comman line; exec_new_project')
        endif
        ! update project info
        call spproj%update_projinfo( cline )
        ! update computer environment
        call spproj%update_compenv( cline )
        ! write project file
        call spproj%write
        ! end gracefully
        call simple_end('**** NEW_PROJECT NORMAL STOP ****')
    end subroutine exec_new_project

    subroutine exec_print_project_info( self, cline )
        class(commander_print_project_info), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        call cline%set('mkdir', 'no')
        call params%new(cline, silent=.true.)
        if(trim(params%json) .eq. 'yes') then
            call spproj%print_info_json(params%projfile) 
        else
            call spproj%print_info(params%projfile)
        endif
        call spproj%kill
    end subroutine exec_print_project_info

    subroutine exec_print_project_vals( self, cline )
        class(commander_print_project_vals), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(binoris)                 :: bos_doc
        type(string)                  :: keys, fname, oritype, str
        logical,          allocatable :: keys_present(:)
        character(len=STDLEN)         :: args(32)
        character(len=STDLEN)         :: str_static
        logical    :: ischar
        type(oris) :: os
        integer    :: nargs, iseg, noris, ikey, iori, state
        real       :: rval, norm(3)
        ! parse the keys
        keys = cline%get_carg('keys')
        str_static = keys%to_char()
        if( str_has_substr(str_static, ',') )then
            call parsestr(str_static, ',', args, nargs)
        else
            args(1) = keys%to_char()
            nargs   = 1
        endif
        ! open the project file
        fname = cline%get_carg('projfile')
        call bos_doc%open(fname)
        ! figure out which segment
        oritype = cline%get_carg('oritype')
        iseg    = oritype2segment(oritype%to_char())
        noris   = bos_doc%get_n_records(iseg)
        if( noris == 0 ) return
        ! read segment
        call os%new(noris, is_ptcl=iseg==3.or.iseg==6) ! see simple_sp_project
        call bos_doc%read_segment(iseg, os)
        ! look for keys
        allocate(keys_present(nargs))
        do ikey=1,nargs
            if( trim(args(ikey)).eq.'eulnorm' )then
                keys_present(ikey) = os%isthere('e1') .and. os%isthere('e2')
            else
                keys_present(ikey) = os%isthere(trim(args(ikey)))
            endif
        end do
        ! print
        if( all(keys_present) )then
            do iori=1,noris
                ! first we write the index
                write(logfhandle,'(i9,a)',advance='no') iori, ' '
                ! then the state
                if( os%isthere(iori,'state') )then
                    state = os%get_state(iori)
                else
                    state = 1
                endif
                write(logfhandle,'(i3,a)',advance='no') state, ' '
                ! then the key values
                do ikey=1,nargs
                    if( trim(args(ikey)).eq.'eulnorm' )then
                        norm = os%get_normal(iori)
                        write(logfhandle,'(f12.4,a)',advance='no') norm(1), ' '
                        write(logfhandle,'(f12.4,a)',advance='no') norm(2), ' '
                        write(logfhandle,'(f12.4,a)',advance='no') norm(3), ' '
                    else
                        ischar = os%ischar(iori,trim(args(ikey)))
                        if( ischar )then
                            call os%getter(iori, trim(args(ikey)), str)
                            write(logfhandle,'(a)',advance='no') str%to_char()//' '
                        else
                            call os%getter(iori, trim(args(ikey)), rval)
                            write(logfhandle,'(f12.4,a)',advance='no') rval, ' '
                        endif
                    endif
                end do
                write(logfhandle,*) ''
            end do
        else
            do ikey=1,nargs
                if( .not. keys_present(ikey) ) write(logfhandle,*) 'key: ', trim(args(ikey)), ' is missing in segment'
            end do
            write(logfhandle,*) 'ERROR! print request failed due to missing keys; simple_commanders_project_core :: exec_print_project_vals'
        endif
    end subroutine exec_print_project_vals

    subroutine exec_print_project_field( self, cline )
        class(commander_print_project_field), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        integer          :: fromto(2)
        logical          :: vol = .false., boxes = .false.
        if(cline%defined("oritype") .and. cline%get_carg('oritype') .eq. 'vol') then
            vol = .true.
            call cline%set("oritype", 'out')
        end if
        call params%new(cline, silent=.true.)
        call spproj%read_segment(params%oritype, params%projfile)
        if(params%boxes .eq. 'yes') boxes = .true.
        if(trim(params%json) .eq. 'yes') then
            if(vol) params%oritype = "vol"
            if(params%fromp .lt. params%top) then
                fromto(1) = params%fromp
                fromto(2) = params%top
                call spproj%print_segment_json(params%oritype, params%projfile, fromto=fromto, sort_key=params%sort, plot_key=params%plot_key, sort_asc=params%sort_asc, hist=params%hist, nran=params%nran, boxes=boxes)
            else
                call spproj%print_segment_json(params%oritype, params%projfile, sort_key=params%sort, plot_key=params%plot_key, sort_asc=params%sort_asc, hist=params%hist, nran=params%nran, boxes=boxes)
            end if
        else
            call spproj%print_segment(params%oritype)
        endif
        call spproj%kill
    end subroutine exec_print_project_field

    subroutine exec_update_project( self, cline )
        class(commander_update_project), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(simple_nice_communicator) :: nice_communicator
        type(parameters)               :: params
        type(sp_project)               :: spproj
        call params%new(cline)
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        ! read relevant segments
        call spproj%read_non_data_segments(params%projfile)
        ! update project info
        call spproj%update_projinfo( cline )
        ! update computer environment
        call spproj%update_compenv( cline )
        ! write the last bit of the project file
        call spproj%write_non_data_segments(params%projfile)
        ! no printing for this program
        call nice_communicator%terminate()
    end subroutine exec_update_project

    subroutine exec_merge_projects( self, cline )
        use simple_commanders_pick, only: commander_extract_distr
        class(commander_merge_projects), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)                :: params
        type(cmdline)                   :: cline_reextract
        type(commander_extract_distr) :: xreextract_distr
        type(sp_project),   allocatable :: spprojs(:)
        integer,            allocatable :: boxes(:), nmics(:)
        real,               allocatable :: smpds(:)
        real    :: smpd
        integer :: nprojs, iproj, box
        logical :: l_reextract, l_has_ptcls
        if( .not.cline%defined('mkdir')        ) call cline%set('mkdir',        'yes')
        if( .not.cline%defined('oritype')      ) call cline%set('oritype',      'ptcl2D')
        if( .not.cline%defined('outside')      ) call cline%set('outside',      'no')
        if( .not.cline%defined('backgr_subtr') ) call cline%set('backgr_subtr', 'no')
        if( .not.cline%defined('pcontrast')    ) call cline%set('pcontrast',    'black')
        ! parameters
        call params%new(cline)
        call cline%set('mkdir','no')
        ! projects info & dimensions
        nprojs = 2 ! only 2 at a time for now
        allocate(spprojs(nprojs),boxes(nprojs),smpds(nprojs),nmics(nprojs))
        call spprojs(1)%read(params%projfile)
        call spprojs(2)%read(params%projfile_target)
        ! all projects have particles?
        l_has_ptcls = spprojs(1)%os_ptcl2D%get_noris() > 0
        do iproj = 2,nprojs
            l_has_ptcls = l_has_ptcls .and. (spprojs(iproj)%os_ptcl2D%get_noris() > 0)
        enddo
        ! gather dimensions
        l_reextract = .false.
        if( l_has_ptcls )then
            do iproj = 1,nprojs
                boxes(iproj) = spprojs(iproj)%get_box()
                smpds(iproj) = spprojs(iproj)%get_smpd()
                nmics(iproj) = spprojs(iproj)%get_nintgs()
            enddo
            l_reextract = any(boxes/=boxes(1)) .or. any(abs(smpds-smpds(1)) > 0.001)
            if( l_reextract .and. any(nmics==0) )then
                THROW_HARD('Cannot merge projects with different particle/pixel sizes without micrographs!')
            endif
        endif
        ! append projects iteratively
        do iproj = 2,nprojs
            call spprojs(1)%append_project(spprojs(iproj))
            call spprojs(iproj)%kill
        enddo
        ! write project
        call spprojs(1)%write(params%projfile)
        call spprojs(1)%kill
        deallocate(spprojs)
        ! reextract/scale particles if necessary
        if( l_reextract )then
            cline_reextract = cline
            call cline_reextract%set('prg', 'reextract')
            if( all(abs(smpds-smpds(1)) < 0.001) )then
                ! all projects have the same pixel size but different particle size
                smpd = smpds(1)
                if( .not.cline%defined('box') )then
                    ! defaults to largest box
                    box = maxval(boxes)
                else
                    box = params%box
                endif
            else
                ! some projects have different pixel size
                ! defaults to largets pixel size
                smpd = maxval(smpds)
                ! defaults to largest particle physical size
                if( .not.cline%defined('box') )then
                    box = floor(maxval(smpds*real(boxes)) / smpd)
                    box = find_larger_magic_box(box)
                else
                    box = params%box
                endif
            endif
            call cline_reextract%set('osmpd', smpd)
            call cline_reextract%set('box',   box)
            ! reextract
            write(logfhandle,'(A,F6.3)')'>>> NEW PIXEL    SIZE: ',smpd
            write(logfhandle,'(A,I6)')  '>>> NEW PARTICLE SIZE: ',box
            call xreextract_distr%execute_safe(cline_reextract)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_PROJECTS NORMAL STOP ****')
    end subroutine exec_merge_projects

    subroutine exec_replace_project_field( self, cline )
        class(commander_replace_project_field), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! projfile <- projfile_target
        call spproj%read(params%projfile)
        call spproj%replace_project(params%projfile_target, params%oritype)
        call spproj%write
        call simple_end('**** REPLACE_PROJECT_FIELD NORMAL STOP ****')
    end subroutine exec_replace_project_field

    subroutine exec_concatenate_projects( self, cline )
        class(commander_concatenate_projects), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        type(string), allocatable :: fnames(:)
        type(sp_project) :: spproj_read, spproj
        type(parameters) :: params   
        integer :: n_spprojs, iproj
        call params%new(cline)
        call read_filetable(params%filetab, fnames)
        n_spprojs = size(fnames)
        call spproj%read(fnames(1))
        do iproj = 2,n_spprojs
            call spproj_read%read(fnames(iproj))
            call spproj%append_project(spproj_read)
        enddo
        call spproj%write(params%projfile)
        call spproj%kill
        call simple_end('**** SIMPLE_CONCATENATE_PROJECTS NORMAL STOP ****')
    end subroutine exec_concatenate_projects

    subroutine exec_aggregate_chunks( self, cline )
        use simple_imgarr_utils,     only: read_cavgs_into_imgarr, dealloc_imgarr
        use simple_strategy2D_utils, only: flag_non_junk_cavgs
        use simple_projfile_utils,   only: merge_chunk_projfiles
        class(commander_aggregate_chunks), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(image),                       allocatable   :: cavg_imgs(:)
        logical,                           allocatable   :: l_non_junk(:)
        real,                              parameter     :: LP_BIN = 20.
        type(sp_project)                                 :: spproj
        type(parameters)                                 :: params
        type(string)                                     :: chunk_fnames(2), folder, cavgs
        real                                             :: smpd, mskrad
        integer                                          :: ldim(3), box, n_non_junk
        ! init params          
        call params%new(cline)
        ! declare constants
        folder          = PATH_HERE
        cavgs           = string('cavgs_merged')//params%ext
        chunk_fnames(1) = params%projfile
        chunk_fnames(2) = params%projfile_target
        ! merge spproj & spproj_target projects
        call merge_chunk_projfiles(chunk_fnames, folder, spproj, write_proj=.false., cavgs_out=cavgs, &
                                   cavgs_replace=.true., sigma2_out=get_fbody(basename(params%projfile), fname2ext(params%projfile)))
        call spproj%write(params%projfile)
        ! read cavgs
        cavg_imgs = read_cavgs_into_imgarr(spproj)
        smpd      = cavg_imgs(1)%get_smpd()
        ldim      = cavg_imgs(1)%get_ldim()
        box       = ldim(1)
        mskrad    = min(real(box/2) - COSMSKHALFWIDTH - 1., 0.5 * params%mskdiam/smpd)
        ! flag non-junk in cavgs
        call flag_non_junk_cavgs( cavg_imgs, LP_BIN, mskrad, l_non_junk )
        ! write verbose_exit_fname file when number non-junk classes >= ncls
        n_non_junk = 0
        if( allocated(l_non_junk) ) n_non_junk = count(l_non_junk)
        if( n_non_junk >= params%ncls .and. params%verbose_exit_fname /= '' ) call simple_touch(params%verbose_exit_fname)
        ! tidy up
        call spproj%kill()
        call dealloc_imgarr(cavg_imgs)
        if( allocated(l_non_junk)   ) deallocate(l_non_junk)
        call simple_end('**** SIMPLE_AGGREGATE_CHUNKS NORMAL STOP ****', verbose_exit=trim(params%verbose_exit).eq.'yes')
    end subroutine exec_aggregate_chunks

    subroutine exec_extract_subproj( self, cline )
        use simple_commanders_project_ptcl, only: commander_import_particles
        class(commander_extract_subproj), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        character(len=*), parameter      :: DEFTAB = 'subproj_deftab.txt'
        integer,          allocatable    :: pinds(:)
        type(string)                     :: ctfflag
        type(parameters)                 :: params
        type(sp_project)                 :: spproj
        type(commander_new_project)      :: xnew_proj
        type(commander_import_particles) :: ximport_particles
        type(ctfparams)                  :: ctfvars
        type(cmdline)                    :: cline_new_proj, cline_import_particles
        type(oris)                       :: os_ptcl2D_prev, os_ptcl3D_prev
        integer :: i, n, np2D, np3D
        call cline%set('mkdir', 'no')
        ! init params
        call params%new(cline)
        ! read the project file
        call spproj%read(params%projfile)
        call spproj%write_segment_inside('projinfo')
        ! create substack
        if( .not. cline%defined('outstk') )then
            params%outstk = params%subprojname//'.mrcs'
        endif
        if( cline%defined('clustind') )then
            call spproj%os_ptcl2D%get_pinds(params%clustind, 'class', pinds)
            n = size(pinds)
            call spproj%write_substk(pinds, params%outstk)
            ! extract previous oris
            os_ptcl2D_prev = spproj%os_ptcl2D%extract_subset(pinds)
            os_ptcl3D_prev = spproj%os_ptcl3D%extract_subset(pinds)
        else
            call spproj%write_substk(params%fromp, params%top, params%outstk)
            ! extract previous oris
            os_ptcl2D_prev = spproj%os_ptcl2D%extract_subset(params%fromp, params%top)
            os_ptcl3D_prev = spproj%os_ptcl3D%extract_subset(params%fromp, params%top)
            ! extract previous pinds
            pinds = nint(os_ptcl3D_prev%get_all('pind'))
            n     = size(pinds) 
        endif
        ! get ctf variables
        ctfvars = spproj%get_ctfparams('stk', 1)
        select case(ctfvars%ctfflag)
            case(CTFFLAG_NO)
                ctfflag = 'no'
            case DEFAULT
                ! generate file with defocus values
                if( cline%defined('class') )then
                    call os_ptcl2D_prev%write(string(DEFTAB))
                else
                    call os_ptcl2D_prev%write(string(DEFTAB))
                endif
                ctfflag = spproj%get_ctfflag('ptcl2D', 1)
        end select
        ! make new project
        call cline_new_proj%set('projname', params%subprojname)
        call xnew_proj%execute_safe(cline_new_proj)
        ! import particles
        call cline_import_particles%set('prg',      'import_particles') ! needs to be here for exec_dir creation
        call cline_import_particles%set('projfile', params%subprojname//'.simple')
        call cline_import_particles%set('cs',       ctfvars%cs)
        call cline_import_particles%set('fraca',    ctfvars%fraca)
        call cline_import_particles%set('kv',       ctfvars%kv)
        call cline_import_particles%set('smpd',     ctfvars%smpd)
        call cline_import_particles%set('stk',      '../'//params%outstk%to_char())
        call cline_import_particles%set('ctf',      ctfflag)
        if( ctfflag .ne. 'no' )then
            call cline_import_particles%set('deftab', '../'//DEFTAB)
        endif
        call ximport_particles%execute_safe(cline_import_particles)
        ! transfer previous particle indices to project
        call spproj%read(params%subprojname//'.simple')
        np3D = spproj%os_ptcl3D%get_noris()
        np2D = spproj%os_ptcl2D%get_noris()
        if( np3D /= n .or. np2D /= n ) THROW_HARD('Incongruent ptcl2D/ptcl3D fields')
        do i = 1,n
            call spproj%os_ptcl2D%transfer_2Dparams(i, os_ptcl2D_prev, i)
            call spproj%os_ptcl3D%transfer_3Dparams(i, os_ptcl3D_prev, i)
            call spproj%os_ptcl2D%set(i, 'pind', pinds(i))
            call spproj%os_ptcl3D%set(i, 'pind', pinds(i))
        end do
        call spproj%write
        ! get back to working dir
        call simple_chdir('../../')
        call del_file(DEFTAB)
        ! destruct
        call spproj%kill
        call os_ptcl2D_prev%kill
        call os_ptcl3D_prev%kill
        ! end gracefully
        call simple_end('**** SIMPLE_EXTRACT_SUBPROJ NORMAL STOP ****')
    end subroutine exec_extract_subproj

    subroutine exec_selection( self, cline )
        class(commander_selection), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters)                :: params
        type(stream_http_communicator)  :: http_communicator
        type(ran_tabu)                  :: rt
        type(sp_project)                :: spproj
        integer,            allocatable :: states(:), ptcls_in_state(:), ptcls_rnd(:)
        integer(kind=kind(ENUM_ORISEG)) :: iseg
        integer                         :: n_lines, fnr, noris, i, nstks, noris_in_state
        integer                         :: state
        logical                         :: l_ctfres, l_icefrac, l_append, l_keep, l_writecls2d, l_writestar
        class(oris), pointer :: pos => NULL()
        l_append     = .false.
        l_writecls2d = .false.
        l_writestar  = .false.
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('prune')  ) call cline%set('prune',  'no')
        if( .not. cline%defined('append') ) call cline%set('append', 'no')
        call params%new(cline, silent=.true.)
        ! http communicator init
        call http_communicator%create(params%niceprocid, params%niceserver%to_char())
        call http_communicator%send_heartbeat()
        if(params%append .eq. 'yes') l_append = .true.
        iseg = oritype2segment(trim(params%oritype))
        ! read project (almost all or largest segments are updated)
        call spproj%read(params%projfile)
        ! update nptcls
        if(.not. cline%defined('nptcls')) params%nptcls = spproj%get_nptcls()
        ! check number of oris in field
        noris = spproj%get_n_insegment(params%oritype)
        ! associate pointer to field
        call spproj%ptr2oritype(params%oritype, pos)
        if( cline%defined('nran') )then
            ! random selection
            if( cline%defined('state') ) then
                call pos%get_pinds(params%state, 'state', ptcls_in_state)
            else
                call pos%get_pinds(1,            'state', ptcls_in_state)
            endif
            noris_in_state = size(ptcls_in_state)
            if( params%nran >= noris_in_state ) THROW_HARD('Random sample size (nran) too small, input a number larger than '//int2str(noris_in_state))
            rt = ran_tabu(noris_in_state)
            allocate(ptcls_rnd(params%nran), source=0)
            call rt%ne_ran_iarr(ptcls_rnd)
            call rt%kill
            ! allocate states and set the state-flags
            allocate(states(noris), source=0)
            do i=1,params%nran
                states(ptcls_rnd(i)) = 1
            enddo
        else if( cline%defined('state') )then
            call pos%get_pinds(params%state, 'state', ptcls_in_state, l_require_updated= .true.)
            noris_in_state = size(ptcls_in_state)
            ! allocate states and set the state-flags
            allocate(states(noris), source=0)
            do i=1,noris_in_state
                states(ptcls_in_state(i)) = 1
            enddo
        else if( cline%defined('ctfresthreshold') .or. cline%defined('icefracthreshold') )then
            l_ctfres  = cline%defined('ctfresthreshold')
            l_icefrac = cline%defined('icefracthreshold')
            if( iseg /= MIC_SEG )then
                THROW_HARD('Incorrect oritype for micrographs selection')
            endif
            allocate(states(noris), source=nint(pos%get_all('state')))
            if( l_ctfres )then
                do i = 1,noris
                    if( states(i) > 0 )then
                        if( pos%isthere(i, 'ctfres') )then
                            if( pos%get(i,'ctfres') > params%ctfresthreshold ) states(i) = 0
                        endif
                    endif
                enddo
            endif
            if( l_icefrac )then
                do i = 1,noris
                    if( states(i) > 0 )then
                        if( pos%isthere(i, 'icefrac') )then
                            if( pos%get(i,'icefrac') > params%icefracthreshold ) states(i) = 0
                        endif
                    endif
                enddo
            endif
        else if( cline%defined('dfmin') )then
            allocate(states(noris), source=nint(pos%get_all('state')))
            do i = 1,noris
                l_keep = .false.
                if( states(i) > 0 )then
                    if( pos%isthere(i, 'dfx') )then
                        if( pos%get(i,'dfx') >= params%dfmin ) l_keep = .true.
                    endif
                    if( pos%isthere(i, 'dfy') )then
                        if( pos%get(i,'dfy') >= params%dfmin ) l_keep = .true.
                    endif
                endif
                if( .not. l_keep ) states(i) = 0
            enddo
        else if( cline%defined('infile') )then
            ! selection based on text file input
            ! sanity check
            n_lines = nlines(params%infile)
            if( cline%defined('state') ) then
                if( spproj%get_n_insegment_state(params%oritype, cline%get_iarg("state")) /= n_lines )then
                    write(logfhandle,*) '# lines in infile '//params%infile%to_char()//': ', n_lines
                    write(logfhandle,*) '# entries in '//trim(params%oritype)//' segment with requested state: ', spproj%get_n_insegment_state(params%oritype, cline%get_iarg("state"))
                    THROW_WARN('# entries in infile/project file '//trim(params%oritype)//' segment with requested state do not match, aborting; exec_selection')
                    return
                endif
            else
                if( noris /= n_lines )then
                    write(logfhandle,*) '# lines in infile '//params%infile%to_char()//': ', n_lines
                    write(logfhandle,*) '# entries in '//trim(params%oritype)//' segment: ', noris
                    THROW_WARN('# entries in infile/project file '//trim(params%oritype)//' segment do not match, aborting; exec_selection')
                    return
                endif
            endif
            ! allocate states and then read the state-flags
            allocate(states(noris))
            call fopen(fnr, FILE=params%infile, STATUS='OLD', action='READ')
            if( cline%defined('state') ) then
                state = cline%get_iarg("state")
                do i=1,noris
                    if( pos%get_state(i) == state ) then
                        read(fnr,*) states(i)
                    else
                        states(i) = 0
                    endif
                end do
            else
                do i=1,n_lines
                    read(fnr,*) states(i)
                end do
            endif
            call fclose(fnr)
        else if( cline%defined('deselfile') )then
            ! selection based on text file containing indices to be deselected
            states=pos%get_all_asint("state")
            n_lines = nlines(params%deselfile)
            call fopen(fnr, FILE=params%deselfile, STATUS='OLD', action='READ')
            do i=1, n_lines
                read(fnr,*) state
                write(logfhandle, *) "DESELECTING INDEX", state
                states(state) = 0
            end do
            call fclose(fnr)
            ! from GUI so write cavgs and star
            l_writecls2d = .true.
            l_writestar  = .true.
        endif
        ! updates relevant segments
        select case(iseg)
            case(MIC_SEG)
                call spproj%report_state2mic(states)
                nstks = spproj%os_stk%get_noris()
                if( nstks > 0 )then
                    call spproj%report_state2stk(states)
                endif
            case(STK_SEG)
                call spproj%report_state2stk(states)
            case(CLS2D_SEG)
                if( trim(params%partition).eq.'yes' )then
                    ! states are partitions of classes
                    call spproj%os_cls2D%set_all('cluster', real(states))
                    where( states /= 0 ) states = 1
                    call spproj%os_cls2D%set_all('state', real(states))
                    ! map partitions to ptcl2D/3D
                    call spproj%map_cls2D_flag_to_ptcls('cluster')
                else
                    call spproj%os_cls2D%set_all('state', real(states))
                    ! map states to ptcl2D/3D & cls3D segments
                    call spproj%map2ptcls_state(append=l_append)
                endif
                if(params%write_imgarr .eq. 'yes') then
                    call spproj%set_cavgs_thumb(params%projfile)
                end if
            case(CLS3D_SEG)
                if(spproj%os_cls3D%get_noris() == spproj%os_cls2D%get_noris())then
                    call spproj%os_cls2D%set_all('state', real(states))
                    ! map states to ptcl2D/3D & cls3D segments
                    call spproj%map2ptcls_state(append=l_append)
                else
                    ! class averages
                    call spproj%os_cls3D%set_all('state', real(states))
                endif
            case(PTCL2D_SEG,PTCL3D_SEG)
                call spproj%os_ptcl2D%set_all('state', real(states))
                call spproj%os_ptcl3D%set_all('state', real(states))
                if( trim(params%prune).eq.'yes' ) call spproj%prune_particles
                call spproj%map_ptcls_state_to_cls
            case DEFAULT
                THROW_HARD('Cannot report selection to segment '//trim(params%oritype)//'; exec_selection')
        end select        
        ! final full write
        call spproj%write(params%projfile)
        if( l_writecls2d ) call spproj%cavgs2mrc()
        if( l_writestar ) then
            if( spproj%os_mic%get_noris() > 0)    call spproj%write_mics_star()
            if( spproj%os_ptcl2D%get_noris() > 0) call spproj%write_ptcl2D_star()
        endif
        call http_communicator%term()
        call simple_end('**** SELECTION NORMAL STOP ****')
    end subroutine exec_selection

end module simple_commanders_project_core
