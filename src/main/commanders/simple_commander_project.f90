! concrete commander: operations on projects (spproject) and associated files
module simple_commander_project
include 'simple_lib.f08'
use simple_binoris_io
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_image,          only: image
use simple_moviewatcher,   only: moviewatcher
use simple_parameters,     only: parameters, params_glob
use simple_sp_project,     only: sp_project
use simple_stack_io,       only: stack_io
use simple_qsys_env,       only: qsys_env
use simple_qsys_funs
implicit none

public :: new_project_commander
public :: print_project_info_commander
public :: print_project_vals_commander
public :: print_project_field_commander
public :: update_project_commander
public :: import_movies_commander
public :: import_boxes_commander
public :: import_particles_commander
public :: import_cavgs_commander
public :: export_cavgs_commander
public :: selection_commander
public :: merge_stream_projects_commander
public :: replace_project_field_commander
public :: scale_project_commander_distr
public :: projops_commander
public :: prune_project_commander_distr
public :: prune_project_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: new_project_commander
  contains
    procedure :: execute      => exec_new_project
end type new_project_commander

type, extends(commander_base) :: print_project_info_commander
  contains
    procedure :: execute      => exec_print_project_info
end type print_project_info_commander

type, extends(commander_base) :: print_project_vals_commander
  contains
    procedure :: execute      => exec_print_project_vals
end type print_project_vals_commander

type, extends(commander_base) :: print_project_field_commander
  contains
    procedure :: execute      => exec_print_project_field
end type print_project_field_commander

type, extends(commander_base) :: update_project_commander
  contains
    procedure :: execute      => exec_update_project
end type update_project_commander

type, extends(commander_base) :: import_movies_commander
  contains
    procedure :: execute      => exec_import_movies
end type import_movies_commander

type, extends(commander_base) :: import_boxes_commander
  contains
    procedure :: execute      => exec_import_boxes
end type import_boxes_commander

type, extends(commander_base) :: import_particles_commander
  contains
    procedure :: execute      => exec_import_particles
end type import_particles_commander

type, extends(commander_base) :: import_cavgs_commander
  contains
    procedure :: execute      => exec_import_cavgs
end type import_cavgs_commander

type, extends(commander_base) :: export_cavgs_commander
  contains
    procedure :: execute      => exec_export_cavgs
end type export_cavgs_commander

type, extends(commander_base) :: selection_commander
  contains
    procedure :: execute      => exec_selection
end type selection_commander

type, extends(commander_base) :: merge_stream_projects_commander
  contains
    procedure :: execute      => exec_merge_stream_projects
end type merge_stream_projects_commander

type, extends(commander_base) :: replace_project_field_commander
  contains
    procedure :: execute      => exec_replace_project_field
end type replace_project_field_commander

type, extends(commander_base) :: scale_project_commander_distr
  contains
    procedure :: execute      => exec_scale_project_distr
end type scale_project_commander_distr

type, extends(commander_base) :: projops_commander
  contains
    procedure :: execute      => exec_projops
end type projops_commander

type, extends(commander_base) :: prune_project_commander_distr
  contains
    procedure :: execute      => exec_prune_project_distr
end type prune_project_commander_distr

type, extends(commander_base) :: prune_project_commander
  contains
    procedure :: execute      => exec_prune_project
end type prune_project_commander

contains

    subroutine exec_new_project( self, cline )
        class(new_project_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        call cline%set('mkdir', 'no')
        call params%new(cline)
        if( cline%defined('projname') .and. cline%defined('dir') )then
            THROW_HARD('both projname and dir defined on command line, use either or; exec_new_project')
        else if( cline%defined('projname') )then
            if( index(params%projname, '.') > 0 ) THROW_HARD('No punctuation allowed in projname')
            if( file_exists(PATH_HERE//trim(params%projname)) )then
                write(logfhandle,*) 'project directory: ', trim(params%projname), ' already exists in cwd: ', trim(params%cwd)
                write(logfhandle,*) 'If you intent to overwrite the existing file, please remove it and re-run new_project'
                THROW_HARD('ABORTING... exec_new_project')
            endif
            ! make project directory
            call simple_mkdir(trim(params%projname), errmsg="commander_project :: new_project;")
            ! change to project directory
            call simple_chdir(trim(params%projname), errmsg="commander_project :: new_project;")
        else if( cline%defined('dir') )then
            params%dir = simple_abspath(trim(params%dir))
            if( .not. file_exists(trim(params%dir)) )then
                write(logfhandle,*) 'input project directory (dir): ', trim(params%dir), ' does not exist'
                THROW_HARD('ABORTING... exec_new_project')
            endif
            call cline%set('projname', basename(params%dir))
            ! change to project directory
            call simple_chdir(trim(params%dir), errmsg="commander_project :: new_project;")
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
        class(print_project_info_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        call cline%set('mkdir', 'no')
        call params%new(cline, silent=.true.)
        call spproj%print_info(params%projfile)
        call spproj%kill
    end subroutine exec_print_project_info

    subroutine exec_print_project_vals( self, cline )
        use simple_sp_project, only: oritype2segment
        class(print_project_vals_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(binoris)                 :: bos_doc
        character(len=:), allocatable :: keys, fname, oritype, str
        logical,          allocatable :: keys_present(:)
        character(len=STDLEN)         :: args(32)
        logical    :: ischar
        type(oris) :: os
        integer    :: nargs, iseg, noris, ikey, iori, state
        real       :: rval, norm(3)
        ! parse the keys
        keys = cline%get_carg('keys')
        if( str_has_substr(keys, ',') )then
            call parsestr(keys, ',', args, nargs)
        else
            args(1) = keys
            nargs   = 1
        endif
        ! open the project file
        fname = cline%get_carg('projfile')
        call bos_doc%open(fname)
        ! figure out which segment
        oritype = cline%get_carg('oritype')
        iseg    = oritype2segment(oritype)
        noris   = bos_doc%get_n_records(iseg)
        if( noris == 0 ) return
        ! read segment
        call os%new(noris, is_ptcl=iseg==3.or.iseg==6) ! see simple_sp_project
        call bos_doc%read_segment(iseg, os)
        ! look for keys
        allocate(keys_present(nargs))
        do ikey=1,nargs
            if( trim(args(ikey)).eq.'eulnorm' )then
                keys_present(ikey) = os%isthere('e1') .and. os%isthere('e1')
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
                    state = nint(os%get(iori,'state'))
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
                            write(logfhandle,'(a)',advance='no') trim(str)//' '
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
            write(logfhandle,*) 'ERROR! print request failed due to missing keys; simple_commander_project :: exec_print_project_vals'
        endif
    end subroutine exec_print_project_vals

    subroutine exec_print_project_field( self, cline )
        class(print_project_field_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        call params%new(cline, silent=.true.)
        call spproj%read_segment(params%oritype, params%projfile)
        call spproj%print_segment(params%oritype)
        call spproj%kill
    end subroutine exec_print_project_field

    subroutine exec_update_project( self, cline )
        class(update_project_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        call params%new(cline)
        ! read relevant segments
        call spproj%read_non_data_segments(trim(params%projfile))
        ! update project info
        call spproj%update_projinfo( cline )
        ! update computer environment
        call spproj%update_compenv( cline )
        ! write the last bit of the project file
        call spproj%write_non_data_segments(trim(params%projfile))
        ! no printing for this program
    end subroutine exec_update_project

    subroutine exec_import_movies( self, cline )
        class(import_movies_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)                       :: params
        type(sp_project)                       :: spproj
        type(oris)                             :: deftab
        type(ctfparams)                        :: ctfvars, prev_ctfvars
        type(moviewatcher)                     :: movie_buff
        character(len=:),          allocatable :: phaseplate
        character(len=LONGSTRLEN), allocatable :: boxfnames(:), movfnames(:)
        character(len=LONGSTRLEN)              :: boxfname
        logical :: inputted_boxtab, inputted_deftab, inputted_dir_movies, inputted_filetab, first_import
        integer :: nmovf, nboxf, i, nprev_movies, nprev_intgs
        ! set defaults
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('ctf')   ) call cline%set('ctf',   'yes')
        call params%new(cline)
        ! parameter input management
        inputted_boxtab = cline%defined('boxtab')
        inputted_deftab = cline%defined('deftab')
        inputted_dir_movies = cline%defined('dir_movies')
        inputted_filetab    = cline%defined('filetab')
        if(inputted_dir_movies .and. inputted_filetab)                        THROW_HARD ('dir_movies cannot be set with a filetab! exec_import_movies')
        if(.not. inputted_dir_movies .and. .not. inputted_filetab)            THROW_HARD ('either dir_movies or filetab must be given! exec_import_movies')
        if(inputted_dir_movies .and. ( inputted_deftab .or. inputted_boxtab)) THROW_HARD ('dir_movies cannot be set with a deftab or boxtab! exec_import_movies')
        ! project file management
        if( .not. file_exists(trim(params%projfile)) )then
            THROW_HARD('project file: '//trim(params%projfile)//' does not exists! exec_import_movies')
        endif
        call spproj%read(params%projfile)
        nprev_intgs  = spproj%get_nintgs()
        nprev_movies = spproj%get_nmovies()
        first_import = nprev_movies==0 .and. nprev_intgs==0
        ! CTF
        if( cline%defined('phaseplate') )then
            phaseplate = cline%get_carg('phaseplate')
        else
            allocate(phaseplate, source='no')
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
            movie_buff = moviewatcher(1)
            call movie_buff%watch( nmovf, movfnames )
            call movie_buff%kill
        endif
        if( params%mkdir.eq.'yes' )then
            ! taking care of paths
            do i=1,nmovf
                if(movfnames(i)(1:1).ne.'/') movfnames(i) = PATH_PARENT//trim(movfnames(i))
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
        if( inputted_boxtab .and. .not. inputted_dir_movies)then
            call read_filetable(params%boxtab, boxfnames)
            nboxf = size(boxfnames)
            if( nboxf /= nmovf )then
                write(logfhandle,*) '# boxfiles: ', nboxf
                write(logfhandle,*) '# movies  : ', nmovf
                THROW_HARD('# boxfiles .ne. # movies; exec_import_movies')
            endif
            do i=1,nmovf
                if( params%mkdir.eq.'yes' )then
                    if(boxfnames(i)(1:1).ne.'/') boxfnames(i) = PATH_PARENT//trim(boxfnames(i))
                endif
                call make_relativepath(CWD_GLOB, boxfnames(i), boxfname)
                if( first_import )then
                    call spproj%os_mic%set_boxfile(i, boxfname)
                else
                    call spproj%os_mic%set_boxfile(nprev_intgs+i, boxfname)
                endif
            end do
        endif 
        ! write project file
        call spproj%write ! full write since projinfo is updated and this is guaranteed to be the first import
        call simple_end('**** IMPORT_MOVIES NORMAL STOP ****')
    end subroutine exec_import_movies

    subroutine exec_import_boxes( self, cline )
        class(import_boxes_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        integer          :: nos_mic, nboxf, i
        character(len=LONGSTRLEN), allocatable :: boxfnames(:)
        character(len=LONGSTRLEN)              :: boxfname
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! project file management
        if( .not. file_exists(trim(params%projfile)) )then
            THROW_HARD('project file: '//trim(params%projfile)//' does not exist! exec_import_boxes')
        endif
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo')
        ! get boxfiles into os_mic
        call read_filetable(params%boxtab, boxfnames)
        nboxf   = size(boxfnames)
        nos_mic = spproj%os_mic%get_noris()
        if( nboxf /= nos_mic )then
            write(logfhandle,*) '# boxfiles       : ', nboxf
            write(logfhandle,*) '# os_mic entries : ', nos_mic
            THROW_HARD('# boxfiles .ne. # os_mic entries; exec_import_boxes')
        endif
        do i=1,nos_mic
            if( params%mkdir.eq.'yes' )then
                if(boxfnames(i)(1:1).ne.'/') boxfnames(i) = PATH_PARENT//trim(boxfnames(i))
            endif
            call make_relativepath(CWD_GLOB,boxfnames(i),boxfname)
            call spproj%os_mic%set_boxfile(i, boxfname)
        end do
        ! write project file
        call spproj%write_segment_inside('mic') ! all that's needed here
        call simple_end('**** IMPORT_BOXES NORMAL STOP ****')
    end subroutine exec_import_boxes

    subroutine exec_import_particles( self, cline )
        class(import_particles_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        character(len=:),      allocatable :: phaseplate, ctfstr
        character(LONGSTRLEN), allocatable :: stkfnames(:)
        real,                  allocatable :: line(:)
        type(parameters) :: params
        type(sp_project) :: spproj
        type(oris)       :: os
        type(nrtxtfile)  :: paramfile
        type(ctfparams)  :: ctfvars
        integer          :: lfoo(3), i, ndatlines, nrecs, n_ori_inputs, nstks
        logical          :: inputted_oritab, inputted_plaintexttab, inputted_deftab
        logical          :: l_stktab_per_stk_parms, is_ptcl
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('ctf')   ) call cline%set('ctf',   'yes')
        l_stktab_per_stk_parms = .true.
        call params%new(cline)
        ! PARAMETER INPUT MANAGEMENT
        ! parameter input flags
        inputted_oritab       = cline%defined('oritab')
        inputted_deftab       = cline%defined('deftab')
        inputted_plaintexttab = cline%defined('plaintexttab')
        n_ori_inputs          = count([inputted_oritab,inputted_deftab,inputted_plaintexttab])
        ! exceptions
        if( n_ori_inputs > 1 )then
            THROW_HARD('multiple parameter sources inputted, please use (oritab|deftab|plaintexttab); exec_import_particles')
        endif
        if( cline%defined('stk') .and. cline%defined('stktab') )then
            THROW_HARD('stk and stktab are both defined on command line, use either or; exec_import_particles')
        endif
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            if( trim(params%ctf) .ne. 'no' )then
                ! there needs to be associated parameters of some form
                if( n_ori_inputs < 1 )then
                    THROW_HARD('stk or stktab input requires associated parameter input when ctf .ne. no (oritab|deftab|plaintexttab)')
                endif
            endif
        else
            THROW_HARD('either stk or stktab needed on command line; exec_import_particles')
        endif
        ! set is particle flag for correct field parsing
        is_ptcl = .false.
        if( cline%defined('stk') ) is_ptcl = .true.
        ! oris input
        if( inputted_oritab )then
            ndatlines = binread_nlines(params%oritab)
            call os%new(ndatlines, is_ptcl=is_ptcl )
            call binread_oritab(params%oritab, spproj, os, [1,ndatlines])
            call spproj%kill ! for safety
        endif
        if( inputted_deftab )then
            ndatlines = binread_nlines(params%deftab)
            call os%new(ndatlines, is_ptcl=is_ptcl )
            call binread_ctfparams_state_eo(params%deftab, spproj, os, [1,ndatlines])
            call spproj%kill ! for safety
        endif
        if( inputted_plaintexttab )then
            call paramfile%new(params%plaintexttab, 1)
            ndatlines = paramfile%get_ndatalines()
            nrecs     = paramfile%get_nrecs_per_line()
            if( nrecs < 1 .or. nrecs > 4 .or. nrecs == 2 )then
                THROW_HARD('unsupported nr of rec:s in plaintexttab; exec_import_particles')
            endif
            call os%new(ndatlines, is_ptcl=is_ptcl )
            allocate( line(nrecs) )
            do i=1,ndatlines
                call paramfile%readNextDataLine(line)
                select case(params%dfunit)
                    case( 'A' )
                        line(1) = line(1)/1.0e4
                        if( nrecs > 1 )  line(2) = line(2)/1.0e4
                    case( 'microns' )
                        ! nothing to do
                    case DEFAULT
                        THROW_HARD('unsupported dfunit; exec_import_particles')
                end select
                select case(params%angastunit)
                    case( 'radians' )
                        if( nrecs == 3 ) line(3) = rad2deg(line(3))
                    case( 'degrees' )
                        ! nothing to do
                    case DEFAULT
                        THROW_HARD('unsupported angastunit; exec_import_particles')
                end select
                select case(params%phshiftunit)
                    case( 'radians' )
                        ! nothing to do
                    case( 'degrees' )
                        if( nrecs == 4 ) line(4) = deg2rad(line(4))
                    case DEFAULT
                        THROW_HARD('unsupported phshiftunit; exec_import_particles')
                end select
                call os%set_dfx(i, line(1))
                if( nrecs > 1 )then
                    call os%set_dfy(i,       line(2))
                    call os%set(i, 'angast', line(3))
                endif
                if( nrecs > 3 )then
                    call os%set(i, 'phshift', line(4))
                endif
            end do
        endif
        if( cline%defined('stktab') )then
            ! importing from stktab
            call read_filetable(params%stktab, stkfnames)
            nstks = size(stkfnames)
            if( params%mkdir.eq.'yes' )then
                do i=1,nstks
                    if(stkfnames(i)(1:1).ne.'/') stkfnames(i) = PATH_PARENT//trim(stkfnames(i))
                    if( .not. file_exists(stkfnames(i)) ) THROW_HARD('modified filetable entry '//trim(stkfnames(i))//' does not exist')
                enddo
            endif
            l_stktab_per_stk_parms = (os%get_noris() == nstks)
            if( (n_ori_inputs == 1) .and. l_stktab_per_stk_parms )then
                ! sampling distance
                call os%set_all2single('smpd', params%smpd)
                ! acceleration voltage
                if( cline%defined('kv') )then
                    call os%set_all2single('kv', params%kv)
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'kv') )then
                            write(logfhandle,*) 'os entry: ', i, ' lacks acceleration volatage (kv)'
                            THROW_HARD('provide kv on command line or update input document; exec_import_particles')
                        endif
                    end do
                endif
                ! spherical aberration
                if( cline%defined('cs') )then
                    call os%set_all2single('cs', params%cs)
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'cs') )then
                            write(logfhandle,*) 'os entry: ', i, ' lacks spherical aberration constant (cs)'
                            THROW_HARD('provide cs on command line or update input document; exec_import_particles')
                        endif
                    end do
                endif
                ! fraction of amplitude contrast
                if( cline%defined('fraca') )then
                    call os%set_all2single('fraca', params%fraca)
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'fraca') )then
                            write(logfhandle,*) 'os entry: ', i, ' lacks fraction of amplitude contrast (fraca)'
                            THROW_HARD('provide fraca on command line or update input document; exec_import_particles')
                        endif
                    end do
                endif
                ! phase-plate
                if( cline%defined('phaseplate') )then
                    call os%set_all2single('phaseplate', trim(params%phaseplate))
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'phaseplate') )then
                            call os%set(i, 'phaseplate', 'no')
                        endif
                    end do
                endif
                call os%getter(1, 'phaseplate', phaseplate)
                if( trim(phaseplate) .eq. 'yes' )then
                    if( .not. os%isthere(1,'phshift') )then
                        THROW_HARD('phaseplate .eq. yes requires phshift input, currently lacking; exec_import_particles')
                    endif
                endif
                ! ctf flag
                if( cline%defined('ctf') )then
                    call os%set_all2single('ctf', trim(params%ctf))
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'ctf') )then
                            call os%set(i, 'ctf', 'yes')
                        endif
                    end do
                endif
                call os%getter(1, 'ctf', ctfstr)
                if( trim(ctfstr) .ne. 'no' )then
                    if( .not. os%isthere(1,'dfx') )then
                        THROW_HARD('ctf .ne. no requires dfx input, currently lacking; exec_import_particles')
                    endif
                endif
            endif
        endif
        if( cline%defined('stk') .or. (cline%defined('stktab').and..not.l_stktab_per_stk_parms) )then
            ctfvars%smpd = params%smpd
            select case(trim(params%ctf))
                case('yes')
                    ctfvars%ctfflag = CTFFLAG_YES
                case('no')
                    ctfvars%ctfflag = CTFFLAG_NO
                case('flip')
                    ctfvars%ctfflag = CTFFLAG_FLIP
                case DEFAULT
                    write(logfhandle,*)
                    THROW_HARD('unsupported ctf flag: '//trim(params%ctf)//'; exec_import_particles')
            end select
            if( ctfvars%ctfflag .ne. CTFFLAG_NO )then
                if( .not. cline%defined('kv')    ) THROW_HARD('kv (acceleration voltage in kV{300}) input required when importing movies; exec_import_particles')
                if( .not. cline%defined('cs')    ) THROW_HARD('cs (spherical aberration constant in mm{2.7}) input required when importing movies; exec_import_particles')
                if( .not. cline%defined('fraca') ) THROW_HARD('fraca (fraction of amplitude contrast{0.1}) input required when importing movies; exec_import_particles')
                if( cline%defined('phaseplate') )then
                    phaseplate = cline%get_carg('phaseplate')
                else
                    allocate(phaseplate, source='no')
                endif
                ctfvars%kv           = params%kv
                ctfvars%cs           = params%cs
                ctfvars%fraca        = params%fraca
                ctfvars%l_phaseplate = phaseplate .eq. 'yes'
            endif
        endif

        ! PROJECT FILE MANAGEMENT
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)

        ! UPDATE FIELDS
        ! add stack if present
        if( cline%defined('stk') )then
            if( n_ori_inputs == 0 .and. trim(params%ctf) .eq. 'no' )then
                ! get number of particles from stack
                call find_ldim_nptcls(params%stk, lfoo, params%nptcls)
                call os%new(params%nptcls, is_ptcl=is_ptcl )
            endif
            ! state = 1 by default
            call os%set_all2single('state', 1.0)
            call os%set_all2single('w',     1.0)
            call spproj%add_single_stk(params%stk, ctfvars, os)
        endif
        ! add list of stacks (stktab) if present
        if( cline%defined('stktab') )then
            ! state = 1 by default
            call os%set_all2single('state', 1.0)
            call os%set_all2single('w',     1.0)
            if( l_stktab_per_stk_parms )then
                ! per stack parameters
                call spproj%add_stktab(stkfnames, os)
            else
                ! per particle parameters
                call spproj%add_stktab(stkfnames, ctfvars, os)
            endif
        endif
        ! WRITE PROJECT FILE
        call spproj%write ! full write since this is guaranteed to be the first import
        call simple_end('**** IMPORT_PARTICLES NORMAL STOP ****')
    end subroutine exec_import_particles

    subroutine exec_import_cavgs( self, cline )
        class(import_cavgs_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        if( file_exists(trim(params%projfile)) ) call spproj%read(params%projfile)
        call spproj%add_cavgs2os_out(params%stk, params%smpd)
        ! update project info
        call spproj%update_projinfo( cline )
        ! update computer environment
        call spproj%update_compenv( cline )
        ! WRITE PROJECT FILE
        call spproj%write ! full write since this is guaranteed to be the first import
        call simple_end('**** IMPORT_CAVGS NORMAL STOP ****')
    end subroutine exec_import_cavgs

    subroutine exec_export_cavgs( self, cline )
        use simple_image, only: image
        class(export_cavgs_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(image)                   :: img
        character(len=:), allocatable :: cavgs_fname, projfile_path
        logical,          allocatable :: lstates(:)
        integer :: ldim(3), icls, ncls, ncavgs, cnt
        real    :: smpd, smpd_phys
        call cline%set('oritype', 'cls2D')
        call params%new(cline)
        ! read files and sanity checks
        if( .not.file_exists(trim(params%projfile)) ) THROW_HARD('Project file does not exist!')
        call spproj%read_segment(params%oritype,params%projfile)
        if( spproj%os_cls2D%get_noris() == 0 ) THROW_HARD('Absent cls2D field!')
        call spproj%read_segment('out',params%projfile)
        lstates = (spproj%os_cls2D%get_all('state') > 0.5)
        if( count(lstates) == 0 ) THROW_HARD('All class averages are deselected')
        call spproj%get_cavgs_stk(cavgs_fname, ncls, smpd)
        if( spproj%os_cls2D%get_noris() /= ncls ) THROW_HARD('Inconsistent # of entries cls2D/out!')
        ! takes care of path
        projfile_path = get_fpath(simple_abspath(params%projfile))
        cavgs_fname   = trim(projfile_path)//'/'//trim(cavgs_fname)
        call find_ldim_nptcls(cavgs_fname, ldim, ncavgs, smpd=smpd_phys)
        if(ncavgs /= ncls)    THROW_HARD('Inconsistent # of cls2D cavgs & physical cavgs!')
        if( abs(smpd-smpd_phys) > 0.001 ) THROW_HARD('Inconsistent sampling distancs in project & physical cavgs!')
        ! copy selected cavgs
        cnt     = 0
        ldim(3) = 1
        call img%new(ldim, smpd)
        do icls = 1,ncls
            if( lstates(icls) )then
                call img%read(cavgs_fname,icls)
                cnt = cnt + 1
                call img%write(params%outstk,cnt)
            endif
        enddo
        ! the end
        call simple_end('**** EXPORT_CAVGS NORMAL STOP ****')
    end subroutine exec_export_cavgs

    subroutine exec_selection( self, cline )
        use simple_sp_project, only: sp_project, oritype2segment
        class(selection_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)                :: params
        type(sp_project)                :: spproj
        integer,            allocatable :: states(:)
        integer(kind=kind(ENUM_ORISEG)) :: iseg
        integer                         :: n_lines,fnr,noris,i,nstks
        real                            :: state
        class(oris), pointer :: pos => NULL()
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('prune') ) call cline%set('prune', 'no')
        call params%new(cline, silent=.true.)
        iseg = oritype2segment(trim(params%oritype))
        ! read project (almost all or largest segments are updated)
        call spproj%read(params%projfile)
        call spproj%update_projinfo( cline )
        ! sanity check
        n_lines = nlines(trim(params%infile))
        noris = spproj%get_n_insegment(params%oritype)
        if( cline%defined('state') ) then
            if( spproj%get_n_insegment_state(params%oritype, cline%get_rarg("state")) /= n_lines )then
                write(logfhandle,*) '# lines in infile '//trim(params%infile)//': ', n_lines
                write(logfhandle,*) '# entries in '//trim(params%oritype)//' segment with requested state: ', noris
                THROW_WARN('# entries in infile/project file '//trim(params%oritype)//' segment with requested state do not match, aborting; exec_selection')
                return
            endif
        else
            noris = spproj%get_n_insegment(params%oritype)
            if( noris /= n_lines )then
                write(logfhandle,*) '# lines in infile '//trim(params%infile)//': ', n_lines
                write(logfhandle,*) '# entries in '//trim(params%oritype)//' segment: ', noris
                THROW_WARN('# entries in infile/project file '//trim(params%oritype)//' segment do not match, aborting; exec_selection')
                return
            endif
        endif
        ! allocate states and then read the state-flags
        allocate(states(noris))
        call fopen(fnr, FILE=trim(params%infile), STATUS='OLD', action='READ')
        if( cline%defined('state') ) then
            state = cline%get_rarg("state")
            call spproj%ptr2oritype(params%oritype, pos)
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
        ! updates relevant segments
        select case(iseg)
            case(MIC_SEG)
                call spproj%os_mic%set_all('state', real(states))
                nstks = spproj%os_stk%get_noris()
                if(nstks > 0)then
                    if( noris /= nstks )then
                        THROW_HARD('This project file has already undergone some selection, use parent project instead')
                    endif
                    call spproj%report_state2stk(states)
                endif
            case(STK_SEG)
                call spproj%report_state2stk(states) ! is this segment update necessary?
            case(CLS2D_SEG)
                call spproj%os_cls2D%set_all('state', real(states))
                call spproj%map2ptcls_state ! map states to ptcl2D/3D & cls3D segments
            case(CLS3D_SEG)
                if(spproj%os_cls3D%get_noris() == spproj%os_cls2D%get_noris())then
                    call spproj%os_cls2D%set_all('state', real(states))
                    call spproj%map2ptcls_state ! map states to ptcl2D/3D & cls3D segments
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
        call simple_end('**** SELECTION NORMAL STOP ****')
    end subroutine exec_selection

    subroutine exec_merge_stream_projects( self, cline )
        class(merge_stream_projects_commander), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj, spproj2merge
        call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! projfile <- projfile_target
        call spproj%read(params%projfile)
        call spproj2merge%read(params%projfile_target)
        call spproj%merge_stream_projects(spproj2merge)
        call spproj%kill
        call simple_end('**** MERGE_STREAMS_PROJECTS NORMAL STOP ****')
    end subroutine exec_merge_stream_projects

    subroutine exec_replace_project_field( self, cline )
        class(replace_project_field_commander), intent(inout) :: self
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

    subroutine exec_scale_project_distr( self, cline )
        use simple_builder,    only: builder
        use simple_parameters, only: params_glob
        class(scale_project_commander_distr), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(chash),      allocatable :: part_params(:)
        character(len=:), allocatable :: projfile_sc
        integer,          allocatable :: parts(:,:)
        integer,          parameter   :: MAX_NCUNITS = 64
        character(len=XLONGSTRLEN)    :: dir_target
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        type(cmdline)    :: cline_scale
        type(parameters) :: params
        type(builder)    :: build
        real             :: smpd, smpd_target
        integer          :: ipart, nparts, nstks, box, newbox, nthr_orig, nparts_orig, ncunits_orig
        logical          :: gen_sc_project
        ! mkdir=yes: a new *_sc project + stacks are generated
        ! mkdir=no : only stacks are scaled
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        gen_sc_project = cline%get_carg('mkdir').eq.'yes'
        ! make parameters and project
        call params%new(cline)
        params%nptcls = 1 ! to avoid excessive memory allocation
        call build%build_spproj(params, cline)
        call build%spproj%read_segment('stk',params%projfile)
        nstks = build%spproj%os_stk%get_noris()
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! copy command line
        cline_scale  = cline
        ! save overridden parameters
        nparts_orig  = params%nparts
        ncunits_orig = params%ncunits
        nthr_orig    = params%nthr
        if( cline%defined('nparts') )then
            nparts = min(params%nparts * params%nthr, nstks)
        else
            nparts = params%nthr
        endif
        call cline_scale%set('nparts', real(nparts))
        smpd = build%spproj%get_smpd()
        box  = build%spproj%get_box()
        if( gen_sc_project )then
            ! make new project & scales
            smpd_target = max(smpd, smpd * real(box)/real(params%newbox))
            dir_target = filepath(PATH_PARENT,'stack_parts_sc')
            if( cline%defined('dir_target') ) dir_target = trim(cline%get_carg('dir_target'))
            call simple_mkdir(dir_target, errmsg="commander_distr_wflows::exec_scale_project_distr ")
            call build%spproj%scale_projfile(smpd_target, projfile_sc, cline, cline_scale, dir=dir_target)
            newbox = nint(cline_scale%get_rarg('newbox'))
            if( newbox == box )then
                write(logfhandle,*)'Inconsistent input dimensions: from ',box,' to ',newbox
                THROW_HARD('inconsistent input dimensions; exec_scale_project_distr')
            endif
        else
            newbox = params%newbox
        endif
        ! needs to be re-set
        call cline_scale%set('smpd', smpd)
        call cline_scale%set('box',  real(box))
        ! setup the environment for distributed execution
        params%nparts       = nparts
        params_glob%nparts  = nparts
        params%ncunits      = min(MAX_NCUNITS, nparts)
        params_glob%ncunits = min(MAX_NCUNITS, nparts)
        params%nthr         = 1
        params_glob%nthr    = 1
        call qenv%new(nparts)
        ! prepares stack-based parts
        parts = split_nobjs_even(nstks, nparts)
        allocate(part_params(nparts))
        do ipart=1,nparts
            call part_params(ipart)%new(2)
            call part_params(ipart)%set('fromp',int2str(parts(ipart,1)))
            call part_params(ipart)%set('top',  int2str(parts(ipart,2)))
        end do
        ! prepare job description
        call cline_scale%gen_job_descr(job_descr)
        call job_descr%set('prg',      'scale')
        call job_descr%set('newbox',   int2str(newbox))
        call job_descr%set('autoscale','no')
        call job_descr%set('nthr',     int2str(1))
        call job_descr%set('nparts',   int2str(nparts))
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr, part_params=part_params, array=L_USE_SLURM_ARR, extra_params=params)
        ! delete copy in working directory
        if( gen_sc_project ) call del_file(params%projfile)
        ! clean
        call qsys_cleanup
        ! end gracefully
        params_glob%nparts  = nparts_orig
        params_glob%ncunits = ncunits_orig
        params_glob%nthr    = nthr_orig
        call build%spproj%kill
        call simple_end('**** SIMPLE_SCALE_PROJECT_DISTR NORMAL STOP ****')
    end subroutine exec_scale_project_distr

    subroutine exec_projops( self, cline )
        class(projops_commander),     intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters)            :: params
        type(sp_project)            :: spproj
        type(oris)                  :: oris_backup
        type(image)                 :: img
        type(stack_io)              :: stksrc, stkdst
        integer, allocatable        :: randmap(:)
        character(len=XLONGSTRLEN)  :: cwd
        integer                     :: i, ifoo, ldim(3)
        logical                     :: l_randomise
        ! init
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call simple_getcwd(cwd)
        call spproj%read( params%projfile )
        l_randomise = trim(params_glob%randomise) .eq. 'yes'
        if( l_randomise ) then
            if(spproj%os_stk%get_noris() .ne. 1)    THROW_HARD('Can only randomise particles in a single stack')
            if(spproj%os_ptcl2D%get_noris() .lt. 1) THROW_HARD('No particle information present in project file')
            write(logfhandle,'(A)')'>>> RANDOMISING PARTICLE ORDER'
            randmap = generate_randomisation_map(spproj%os_ptcl2D%get_noris(), 5)
            write(logfhandle,'(A)')'>>> REMAPPING PARTICLES'
            oris_backup = spproj%os_ptcl2D
            do i=1, spproj%os_ptcl2D%get_noris()
                call spproj%os_ptcl2D%transfer_ori(i, oris_backup, randmap(i))
            enddo
            call oris_backup%kill
            write(logfhandle,'(A)')'>>> WRITING UPDATED STACK'
            if(.not. file_exists(spproj%os_stk%get_static(1, 'stk'))) THROW_HARD('Stack file does not exist')
            call find_ldim_nptcls(spproj%os_stk%get_static(1, 'stk'), ldim, ifoo)
            ldim(3) = 1
            call img%new(ldim, params%smpd)
            call stkdst%open(trim(params%outstk), params%smpd, 'write', box=ldim(1))
            call stksrc%open(spproj%os_stk%get_static(1, 'stk'), params%smpd, 'read', bufsz=spproj%os_ptcl2D%get_noris())
            call stksrc%read_whole ! can we make this better/faster/less ram intensive?
            do i=1, spproj%os_ptcl2D%get_noris()
                call stksrc%read(randmap(i), img)
                call stkdst%write(i, img)
            enddo
            call stkdst%close
            call stksrc%close
            call img%kill
            spproj%os_ptcl2D = spproj%os_ptcl3D
            call spproj%os_stk%set(1,'stk', trim(cwd) // '/' // trim(params%outstk))
        endif
        ! update project info
        call spproj%update_projinfo( cline )
        ! update computer environment
        call spproj%update_compenv( cline )
        ! write project file
        call spproj%write(basename(params%projfile))
        ! end gracefully
        if (allocated(randmap)) deallocate(randmap)
        call simple_end('**** PROJOPS NORMAL STOP ****')
        
        contains
        
            function generate_randomisation_map( array_size, niter ) result( array )
                integer, intent(in)    :: array_size, niter
                integer, allocatable   :: array(:)
                integer                :: i, iswap, iter, tmp
                real                   :: rrand
                array = [( i, i=1, array_size )]
                do iter=1,niter
                    do i=1, array_size
                        call random_number(rrand)
                        iswap = floor( array_size * rrand) + 1
                        tmp   = array(iswap)
                        array(iswap) = array(i)
                        array(i)     = tmp
                    enddo
                enddo
                
            end function generate_randomisation_map
            
        
    end subroutine exec_projops

    subroutine exec_prune_project_distr( self, cline )
        class(prune_project_commander_distr), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(cmdline)                 :: cline_distr
        type(sp_project)              :: spproj
        type(sp_project), allocatable :: spproj_part(:)
        type(qsys_env)                :: qenv
        type(chash)                   :: job_descr
        type(chash),      allocatable :: part_params(:)
        integer,          allocatable :: parts(:,:), pinds(:)
        real,             allocatable :: states(:)
        character(len=:), allocatable :: fname
        logical,          allocatable :: part_mask(:)
        integer :: imic,nmics,cnt,istk,nstks,ipart,nptcls,nparts,iptcl
        integer :: nstks_orig,nptcls_orig,nmics_orig,i
        ! init
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! sanity checks
        call spproj%read(params%projfile)
        nstks = spproj%get_n_insegment('stk')
        if( nstks == 0 ) THROW_HARD('No stack to process!')
        nptcls = spproj%get_n_insegment('ptcl2D')
        if( nstks == 0 ) THROW_HARD('No particles to process!')
        ! identify the particle indices with state .ne. 0 unless state is given on command line
        states = spproj%os_ptcl2D%get_all('state')
        allocate(pinds(nptcls), source=(/(i,i=1,nptcls)/))
        if( cline%defined('state') )then
            pinds = pack(pinds, mask=abs(states - real(params%state)) < 0.1)
        else
            pinds = pack(pinds, mask=states > 0.5)
        endif
        call arr2txtfile(pinds, 'selected_indices.txt')
        deallocate(states, pinds)
        nmics       = spproj%get_nintgs()
        nstks_orig  = nstks
        nptcls_orig = nptcls
        nmics_orig  = nmics
        ! DISTRIBUTED EXECUTION
        ! setup the environment for distributed execution
        cline_distr = cline
        call cline_distr%set('prg',  'prune_project')
        call cline_distr%set('nthr', 1.)
        call cline_distr%set('oritype', 'stk')
        nparts = min(params%nparts, nstks)
        allocate(spproj_part(nparts),part_mask(nparts))
        parts     = split_nobjs_even(nstks, nparts)
        part_mask = .true.
        allocate(part_params(nparts))
        do ipart=1,nparts
            call part_params(ipart)%new(2)
            call part_params(ipart)%set('fromp',int2str(parts(ipart,1)))
            call part_params(ipart)%set('top',  int2str(parts(ipart,2)))
        end do
        call qenv%new(nparts)
        ! prepare job description
        call cline_distr%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(job_descr, part_params=part_params, array=L_USE_SLURM_ARR)
        ! ASSEMBLY
        do ipart = 1,nparts
            fname = trim(ALGN_FBODY)//int2str(ipart)//METADATA_EXT
            part_mask(ipart) = file_exists(fname)
        enddo
        ! copy updated micrographs
        nstks = 0
        if( nmics > 0 )then
            cnt = 0
            do ipart = 1,nparts
                if( .not.part_mask(ipart) ) cycle
                fname = trim(ALGN_FBODY)//int2str(ipart)//METADATA_EXT
                call spproj_part(ipart)%read_segment('mic',fname)
                cnt = cnt + spproj_part(ipart)%os_mic%get_noris()
            enddo
            nmics = cnt
            if( nmics > 0 )then
                call spproj%os_mic%new(nmics, is_ptcl=.false.)
                cnt = 0
                do ipart = 1,nparts
                    if( spproj_part(ipart)%os_mic%get_noris() == 0 ) cycle
                    do imic = 1,spproj_part(ipart)%os_mic%get_noris()
                        cnt = cnt + 1
                        call spproj%os_mic%transfer_ori(cnt, spproj_part(ipart)%os_mic, imic)
                        if( nint(spproj%os_mic%get(cnt,'nptcls')) > 0 ) nstks = nstks + 1
                    enddo
                    call spproj_part(ipart)%kill
                enddo
                if( nstks /= nmics ) THROW_HARD('Inconsistent number of stacks and micrographs!')
            endif
            write(logfhandle,'(A,I8,A,I8)')'>>> # OF MICROGRAPHS:',nmics_orig,' -> ', nmics
        endif
        ! copy updated stacks
        cnt = 0
        do ipart = 1,nparts
            if( .not.part_mask(ipart) ) cycle
            fname = trim(ALGN_FBODY)//int2str(ipart)//METADATA_EXT
            call spproj_part(ipart)%read_segment('stk',fname)
            cnt = cnt + spproj_part(ipart)%os_stk%get_noris()
        enddo
        nstks = cnt
        write(logfhandle,'(A,I8,A,I8)')'>>> # OF STACKS     :',nstks_orig,' -> ', nstks
        if( nstks > 0 )then
            call spproj%os_stk%new(nstks, is_ptcl=.false.)
            cnt = 0
            do ipart = 1,nparts
                if( spproj_part(ipart)%os_stk%get_noris() == 0 ) cycle
                do istk = 1,spproj_part(ipart)%os_stk%get_noris()
                    cnt = cnt + 1
                    call spproj%os_stk%transfer_ori(cnt, spproj_part(ipart)%os_stk, istk)
                enddo
                call spproj_part(ipart)%kill
            enddo
            ! copy updated particles 2D segment
            cnt = 0
            do ipart = 1,nparts
                if( .not.part_mask(ipart) ) cycle
                fname = trim(ALGN_FBODY)//int2str(ipart)//METADATA_EXT
                call spproj_part(ipart)%read_segment('ptcl2D',fname)
                cnt = cnt + spproj_part(ipart)%os_ptcl2D%get_noris()
            enddo
            nptcls = cnt
            call spproj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
            cnt = 0
            do ipart = 1,nparts
                if( spproj_part(ipart)%os_ptcl2D%get_noris() == 0 ) cycle
                do iptcl = 1,spproj_part(ipart)%os_ptcl2D%get_noris()
                    cnt = cnt + 1
                    call spproj%os_ptcl2D%transfer_ori(cnt, spproj_part(ipart)%os_ptcl2D, iptcl)
                enddo
                call spproj_part(ipart)%kill
            enddo
            ! copy updated particles 3D segment
            call spproj%os_ptcl3D%new(nptcls, is_ptcl=.true.)
            cnt = 0
            do ipart = 1,nparts
                if( .not.part_mask(ipart) ) cycle
                fname = trim(ALGN_FBODY)//int2str(ipart)//METADATA_EXT
                call spproj_part(ipart)%read_segment('ptcl3D',fname)
                if( spproj_part(ipart)%os_ptcl3D%get_noris() == 0 ) cycle
                do iptcl = 1,spproj_part(ipart)%os_ptcl3D%get_noris()
                    cnt = cnt + 1
                    call spproj%os_ptcl3D%transfer_ori(cnt, spproj_part(ipart)%os_ptcl3D, iptcl)
                enddo
                call spproj_part(ipart)%kill
            enddo
        endif
        write(logfhandle,'(A,I8,A,I8)')'>>> # OF PARTICLES  :',nptcls_orig,' -> ', nptcls
        ! final write
        call spproj%write(params%projfile)
        ! clean up
        call spproj%kill
        call qsys_cleanup
        do ipart = 1,nparts
            call part_params(ipart)%kill
            call spproj_part(ipart)%kill
            call del_file(trim(ALGN_FBODY)//int2str(ipart)//METADATA_EXT)
        enddo
        deallocate(spproj_part,part_params)
        ! end gracefully
        call simple_end('**** SIMPLE_PRUNE_PROJECT_DISTR NORMAL STOP ****')
    end subroutine exec_prune_project_distr

    subroutine exec_prune_project( self, cline )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_qsys_funs, only: qsys_job_finished
        use simple_image,     only: image
        use simple_cmdline,   only: cmdline
        class(prune_project_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)              :: params
        type(image)                   :: img
        type(sp_project)              :: spproj, spproj_out
        type(ori)                     :: o_stk
        character(len=:), allocatable :: newstkname,stkname,ext
        logical,          allocatable :: stks_mask(:), ptcls_mask(:)
        integer,          allocatable :: stkinds(:), stk2mic_inds(:), mic2stk_inds(:)
        character(len=LONGSTRLEN) :: relstkname, stkdir
        real                      :: smpd
        integer                   :: iptcl, istk, stk_cnt, nptcls_tot, ptcl_cnt
        integer                   :: box, nstks, nstks_tot, fromp, top, fromp_glob, top_glob, nmics_tot
        integer                   :: nstks_part, nptcls_part, stkind, nstks_prev, ptcl_glob
        ! init
        call params%new(cline)
        if( params%dir .eq. '' )then
            stkdir = PATH_HERE
        else
            stkdir = trim(params%dir)
        endif
        ! particles
        call spproj%read_segment('ptcl2D', params%projfile)
        nptcls_tot = spproj%os_ptcl2D%get_noris()
        allocate(ptcls_mask(nptcls_tot), stkinds(nptcls_tot))
        nptcls_part = 0
        !$omp parallel do proc_bind(close) default(shared) private(iptcl) reduction(+:nptcls_part)
        do iptcl=1,nptcls_tot
            ptcls_mask(iptcl) = spproj%os_ptcl2D%get_state(iptcl) > 0
            if( ptcls_mask(iptcl) )then
                stkinds(iptcl) = nint(spproj%os_ptcl2D%get(iptcl,'stkind'))
                if( stkinds(iptcl) >= params%fromp .and. stkinds(iptcl) <= params%top )then
                    nptcls_part = nptcls_part+1
                endif
            else
                stkinds(iptcl) = 0
            endif
        enddo
        !$omp end parallel do
        call spproj%read_segment('ptcl3D', params%projfile)
        call spproj_out%os_ptcl2D%new(nptcls_part, is_ptcl=.true.)
        call spproj_out%os_ptcl3D%new(nptcls_part, is_ptcl=.true.)
        ! stacks
        call spproj%read_segment('stk', params%projfile)
        nstks_tot = spproj%get_nstks()
        if( nstks_tot == 0 ) THROW_HARD('No images to operate on!')
        allocate(stks_mask(nstks_tot))
        do istk=1,nstks_tot
            stks_mask(istk) = spproj%os_stk%get_state(istk) > 0
            if( count(stkinds==istk) == 0 ) stks_mask(istk) = .false.
        enddo
        nstks = count(stks_mask)
        nstks_part = count(stks_mask(params%fromp:params%top))
        if( nstks_part == 0 )then
            call qsys_job_finished('simple_commander_project :: exec_prune_project')
            return
        endif
        call spproj_out%os_stk%new(nstks_part, is_ptcl=.false.)
        ! micrographs
        call spproj%read_segment('mic', params%projfile)
        nmics_tot = spproj%os_mic%get_noris()
        if( nmics_tot > 0 )then
            call spproj%get_mic2stk_inds(mic2stk_inds, stk2mic_inds)
            call spproj_out%os_mic%new(nstks_part, is_ptcl=.false.)
        endif
        ! new stacks
        box  = spproj%get_box()
        smpd = spproj%get_smpd()
        write(logfhandle,'(A)')'>>> GENERATING STACK(S)'
        call img%new([box,box,1],smpd)
        call simple_mkdir(stkdir)
        nstks_prev = count(stks_mask(:params%fromp-1))
        stkind     = nstks_prev
        stk_cnt    = 0
        if( params%fromp == 1 )then
            top_glob   = 0
        else
            top = nint(spproj%os_stk%get(params%fromp-1,'top'))
            top_glob  = count(ptcls_mask(1:top))
        endif
        ptcl_glob  = 0
        do istk=params%fromp,params%top
            if( .not.stks_mask(istk) ) cycle
            stk_cnt = stk_cnt + 1
            stkind  = stkind  + 1
            call spproj%os_stk%get_ori(istk, o_stk)
            call o_stk%getter('stk',stkname)
            ext        = fname2ext(stkname)
            newstkname = trim(stkdir)//trim(get_fbody(basename(stkname),ext))//trim(STK_EXT)
            fromp      = nint(o_stk%get('fromp'))
            top        = nint(o_stk%get('top'))
            fromp_glob = top_glob+1
            ptcl_cnt   = 0
            do iptcl=fromp,top
                if( .not.ptcls_mask(iptcl) )cycle
                ptcl_glob = ptcl_glob + 1
                top_glob  = top_glob+1
                ptcl_cnt  = ptcl_cnt+1
                ! copy image
                if(spproj%os_ptcl2D%isthere(iptcl, 'indstk') .and. spproj%os_ptcl2D%get(iptcl, 'indstk') > 0.0) then
                        write(logfhandle, *) "STK ", spproj%os_ptcl2D%get(iptcl,'indstk')
                        call img%read(stkname, nint(spproj%os_ptcl2D%get(iptcl,'indstk')))
                else
                        write(logfhandle, *) "STK2 " // int2str(nint(spproj%os_ptcl2D%get(iptcl,'indstk')))
                        call img%read(stkname, iptcl-fromp+1)
                endif
                call img%write(newstkname, ptcl_cnt)
                ! update orientations
                call spproj_out%os_ptcl2D%transfer_ori(ptcl_glob, spproj%os_ptcl2D, iptcl)
                call spproj_out%os_ptcl3D%transfer_ori(ptcl_glob, spproj%os_ptcl3D, iptcl)
                call spproj_out%os_ptcl2D%set(ptcl_glob,'stkind',real(stkind))
                call spproj_out%os_ptcl3D%set(ptcl_glob,'stkind',real(stkind))
                call spproj_out%os_ptcl2D%set(ptcl_glob,'indstk',real(ptcl_cnt))
                call spproj_out%os_ptcl3D%set(ptcl_glob,'indstk',real(ptcl_cnt))
            enddo
            ! update stack
            call make_relativepath(CWD_GLOB, newstkname, relstkname)
            call o_stk%set('stk',   relstkname)
            call o_stk%set('fromp', real(fromp_glob))
            call o_stk%set('top',   real(top_glob))
            call o_stk%set('nptcls',real(ptcl_cnt))
            call spproj_out%os_stk%set_ori(stk_cnt, o_stk)
            ! update micrograph
            if( nmics_tot > 0 ) then
                call spproj_out%os_mic%transfer_ori(stk_cnt, spproj%os_mic, stk2mic_inds(istk))
                call spproj_out%os_mic%set(stk_cnt,'nptcls',real(ptcl_cnt))
            endif
        enddo
        spproj_out%projinfo = spproj%projinfo
        spproj_out%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris() > 0 ) spproj_out%jobproc = spproj%jobproc
        call spproj%kill
        call spproj_out%write(trim(ALGN_FBODY)//int2str(params%part)//METADATA_EXT)
        ! cleanup
        call spproj_out%kill
        call img%kill
        call o_stk%kill
        ! end gracefully
        call qsys_job_finished('simple_commander_project :: exec_prune_project')
    end subroutine exec_prune_project

end module simple_commander_project
