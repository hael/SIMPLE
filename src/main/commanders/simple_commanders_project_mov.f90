!@descr: project commanders for movie-related things
module simple_commanders_project_mov
use simple_commander_module_api
use simple_stream_watcher, only: stream_watcher
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_import_movies
  contains
    procedure :: execute      => exec_import_movies
end type commander_import_movies

type, extends(commander_base) :: commander_write_mic_filetab
  contains
    procedure :: execute      => exec_write_mic_filetab
end type commander_write_mic_filetab

contains

    subroutine exec_import_movies( self, cline )
        class(commander_import_movies), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(simple_nice_communicator) :: nice_communicator
        type(parameters)               :: params
        type(sp_project)               :: spproj
        type(oris)                     :: deftab
        type(ctfparams)                :: ctfvars, prev_ctfvars
        type(stream_watcher)           :: movie_buff
        type(string)                   :: phaseplate, boxfname
        type(string), allocatable      :: boxfnames(:), movfnames(:) 
        logical :: inputted_boxtab, inputted_deftab, inputted_dir_movies, inputted_filetab, first_import
        integer :: nmovf, nboxf, i, nprev_movies, nprev_intgs
        ! set defaults
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('ctf')   ) call cline%set('ctf',   'yes')
        call params%new(cline)
        ! parameter input management
        inputted_boxtab     = cline%defined('boxtab')
        inputted_deftab     = cline%defined('deftab')
        inputted_dir_movies = cline%defined('dir_movies')
        inputted_filetab    = cline%defined('filetab')
        if(inputted_dir_movies .and. inputted_filetab)                        THROW_HARD ('dir_movies cannot be set with a filetab! exec_import_movies')
        if(.not. inputted_dir_movies .and. .not. inputted_filetab)            THROW_HARD ('either dir_movies or filetab must be given! exec_import_movies')
        if(inputted_dir_movies .and. ( inputted_deftab .or. inputted_boxtab)) THROW_HARD ('dir_movies cannot be set with a deftab or boxtab! exec_import_movies')
        ! project file management
        if(.not. file_exists(params%projfile))then
            THROW_HARD('project file: '//params%projfile%to_char()//' does not exists! exec_import_movies')
        endif
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        call spproj%read(params%projfile)
        nprev_intgs  = spproj%get_nintgs()
        nprev_movies = spproj%get_nmovies()
        first_import = nprev_movies==0 .and. nprev_intgs==0
        ! CTF
        if( cline%defined('phaseplate') )then
            phaseplate = cline%get_carg('phaseplate')
        else
            phaseplate ='no'
        endif
        ctfvars%smpd  = params%smpd
        ctfvars%kv    = params%kv
        ctfvars%cs    = params%cs
        ctfvars%fraca = params%fraca
        select case(params%ctf)
            case('yes')
                ctfvars%ctfflag = CTFFLAG_YES
            case('no')
                ctfvars%ctfflag = CTFFLAG_NO
            case('flip')
                ctfvars%ctfflag = CTFFLAG_FLIP
            case DEFAULT
                THROW_HARD('ctf flag: '//trim(params%ctf)//' not supported; exec_import_movies')
        end select
        ctfvars%l_phaseplate = .false.
        if( trim(params%phaseplate) .eq. 'yes' ) ctfvars%l_phaseplate = .true.
        ! sanity checks
        if( first_import )then
            if( spproj%get_nstks()>0)then
                THROW_HARD('Can only import movies/micrographs to an empty project! exec_import_movies')
            endif
        else
            if( spproj%get_nstks()>0)then
                THROW_HARD('Can only import movies/micrographs to a project with movies/micrographs only! exec_import_movies')
            endif
            if( inputted_boxtab )then
                if( .not.spproj%has_boxfile() )then
                    THROW_HARD('Can only import boxes to a project with existing boxes! exec_import_movies')
                endif
            endif
            prev_ctfvars = spproj%get_micparams(1)
            if( ctfvars%ctfflag /= prev_ctfvars%ctfflag ) THROW_HARD('CTF infos do not match! exec_import_movies')
            if( .not.is_equal(ctfvars%smpd, prev_ctfvars%smpd)  ) THROW_HARD('The sampling distances do not match! exec_import_movies')
            if( .not.is_equal(ctfvars%cs,   prev_ctfvars%cs)   ) THROW_HARD('The spherical aberrations do not match! exec_import_movies')
            if( .not.is_equal(ctfvars%kv,   prev_ctfvars%kv)   ) THROW_HARD('The voltages do not match! exec_import_movies')
            if( .not.is_equal(ctfvars%fraca,prev_ctfvars%fraca)) THROW_HARD('The amplitude contrasts do not match! exec_import_movies')
            if( ctfvars%l_phaseplate.neqv.prev_ctfvars%l_phaseplate ) THROW_HARD('Phaseplate infos do not match! exec_import_movies')
        endif
        ! update project info
        call spproj%update_projinfo( cline )
        ! updates segment
        nmovf = 0
        if( inputted_filetab ) then
            nmovf = nlines(params%filetab)
            call read_filetable(params%filetab, movfnames)
        else if( inputted_dir_movies) then
            ! movie watcher init for files older that 1 second (assumed already in place at exec)
            movie_buff = stream_watcher(1,params%dir_movies)
            call movie_buff%watch( nmovf, movfnames )
            call movie_buff%kill
        endif
        if( params%mkdir.eq.'yes' )then
            ! taking care of paths
            do i=1,nmovf
                if(movfnames(i)%to_char([1,1]).ne.'/') movfnames(i) = PATH_PARENT//movfnames(i)%to_char()
            enddo
        endif
        if( inputted_deftab .and. .not. inputted_dir_movies )then
            ! micrographs with pre-determined CTF parameters
            call deftab%new(nlines(params%deftab), is_ptcl=.false.)
            call deftab%read_ctfparams_state_eo(params%deftab)
            call spproj%add_intgs(movfnames, deftab, ctfvars)
            call deftab%kill
        else
            ! movies/micrographs
            call spproj%add_movies(movfnames, ctfvars)
        endif
        ! add boxtab
        if( inputted_boxtab .and. .not. inputted_dir_movies )then
            call read_filetable(params%boxtab, boxfnames)
            nboxf = size(boxfnames)
            if( nboxf /= nmovf )then
                write(logfhandle,*) '# boxfiles: ', nboxf
                write(logfhandle,*) '# movies  : ', nmovf
                THROW_HARD('# boxfiles .ne. # movies; exec_import_movies')
            endif
            do i=1,nmovf
                if( trim(params%mkdir).eq.'yes' )then
                    if(boxfnames(i)%to_char([1,1]).ne.'/') boxfnames(i) = PATH_PARENT//boxfnames(i)%to_char()
                endif
                boxfname = simple_abspath(boxfnames(i))
                if( first_import )then
                    call spproj%set_boxfile(i, boxfname)
                else
                    call spproj%set_boxfile(nprev_intgs+i, boxfname)
                endif
            end do
        endif 
        ! write project file
        call spproj%write ! full write since projinfo is updated and this is guaranteed to be the first import
        call nice_communicator%terminate()
        call simple_end('**** IMPORT_MOVIES NORMAL STOP ****')
    end subroutine exec_import_movies

    subroutine exec_write_mic_filetab( self, cline )
        class(commander_write_mic_filetab), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(string), allocatable :: micstab(:)
        integer,      allocatable :: orimap(:)
        type(parameters) :: params
        type(sp_project) :: spproj
        if( .not.cline%defined('mkdir') ) call cline%set('mkdir', 'no')
        call params%new(cline)
        call spproj%read(params%projfile)
        call spproj%get_mics_table(micstab, orimap)
        call write_filetable(params%fname, micstab)
        call spproj%kill
        if( allocated(micstab) ) deallocate(micstab)
        if( allocated(orimap)  ) deallocate(orimap)
        call simple_end('**** WRITE_MIC_FILETAB NORMAL STOP ****')
    end subroutine exec_write_mic_filetab

end module simple_commanders_project_mov
