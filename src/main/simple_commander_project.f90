! concrete commander: operations on projects (spproject) and associated files
module simple_commander_project
include 'simple_lib.f08'
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_sp_project,     only: sp_project
use simple_oris,           only: oris
use simple_binoris_io,     only: binread_nlines, binread_oritab
use simple_parameters,     only: parameters
implicit none

public :: project2txt_commander
public :: txt2project_commander
public :: print_project_info_commander
public :: report_selection_commander
public :: print_project_vals_commander
public :: print_project_field_commander
public :: new_project_commander
public :: update_project_commander
public :: import_movies_commander
public :: import_boxes_commander
public :: import_particles_commander
public :: import_cavgs_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: project2txt_commander
  contains
    procedure :: execute      => exec_project2txt
end type project2txt_commander
type, extends(commander_base) :: txt2project_commander
  contains
    procedure :: execute      => exec_txt2project
end type txt2project_commander
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
type, extends(commander_base) :: report_selection_commander
  contains
    procedure :: execute      => exec_report_selection
end type report_selection_commander
type, extends(commander_base) :: new_project_commander
  contains
    procedure :: execute      => exec_new_project
end type new_project_commander
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

contains

    !> convert text (.txt) oris doc to binary (.simple)
    subroutine exec_txt2project( self, cline )
        class(txt2project_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters) :: params
        type(oris)       :: os
        type(sp_project) :: spproj
        integer          :: noris
        logical          :: projfile_existed
        call params%new(cline)
        noris = nlines(params%oritab)
        call os%new(noris)
        call os%read(params%oritab)
        projfile_existed = .false.
        if( file_exists(params%projfile) )then
            call spproj%read(params%projfile)
            projfile_existed = .true.
        endif
        call spproj%set_sp_oris(params%oritype, os)
        if( projfile_existed )then
            call spproj%write_segment_inside(params%oritype, params%projfile)
        else
            call spproj%write(params%projfile)
        endif
        call spproj%kill
        call simple_end('**** TXT2PROJECT NORMAL STOP ****')
    end subroutine exec_txt2project

    !> convert binary (.simple) oris doc to text (.txt)
    subroutine exec_project2txt( self, cline )
        class(project2txt_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        call params%new(cline)
        call spproj%read_segment(params%oritype, params%projfile)
        call spproj%write_segment2txt(params%oritype, params%outfile)
        call spproj%kill
        call simple_end('**** PROJECT2TXT NORMAL STOP ****')
    end subroutine exec_project2txt

    !> convert binary (.simple) oris doc to text (.txt)
    subroutine exec_print_project_info( self, cline )
        class(print_project_info_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        call params%new(cline, silent=.true.)
        call spproj%print_info(params%projfile)
        call spproj%kill
        ! no additional printing
    end subroutine exec_print_project_info

    !> prints the values of inputted keys in the inputted segment
    subroutine exec_print_project_vals( self, cline )
        use simple_binoris,    only: binoris
        use simple_sp_project, only: oritype2segment
        use simple_binoris,    only: binoris
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
        call os%new(noris)
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

    subroutine exec_report_selection( self, cline )
        use simple_binoris,    only: binoris
        use simple_sp_project, only: sp_project, oritype2segment
        use simple_oris,       only: oris
        class(report_selection_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(binoris)    :: bos_doc
        type(oris)       :: os
        type(parameters) :: params
        type(sp_project) :: spproj
        integer, allocatable :: states(:)
        integer(kind=kind(ENUM_ORISEG)) :: iseg
        integer :: n_lines, fnr, noris, i
        call params%new(cline, silent=.true.)
        ! read the state-flags
        n_lines = nlines(trim(params%infile))
        allocate(states(n_lines))
        call fopen(fnr, FILE=trim(params%infile), STATUS='OLD', action='READ')
        do i=1,n_lines
            read(fnr,*) states(i)
        end do
        call fclose(fnr)
        iseg = oritype2segment(trim(params%oritype))
        if( iseg == CLS2D_SEG )then
            ! different approach since state mapping may have to be propagated to ptcl fields
            call spproj%read(params%projfile) ! full read since multiple segments will be modified
            noris = spproj%os_cls2D%get_noris()
            if( noris /= n_lines )then
                write(logfhandle,*) '# lines in infile         : ', n_lines
                write(logfhandle,*) '# entries in CLS2D segment: ', noris
                THROW_WARN('# entries in infile/project file CLS2D segment do not match, aborting; exec_report_selection')
                return
            endif
            ! update states
            call spproj%os_cls2D%set_all('state', real(states))
            ! map states to ptcl segments
            call spproj%map2ptcls_state
            call spproj%write(params%projfile) ! full write since multiple segments possibly modified
        else
            ! look in projfile
            call bos_doc%open(trim(params%projfile))
            noris = bos_doc%get_n_records(iseg)
            if( noris == 0 )then
                call bos_doc%close
                THROW_WARN('empty project file segment, nothing to update, aborting; exec_report_selection')
                return
            endif
            if( noris /= n_lines )then
                call bos_doc%close
                write(logfhandle,*) '# lines in infile        : ', n_lines
                write(logfhandle,*) '# entries in file segment: ', noris
                THROW_WARN('# entries in infile/project file segment do not match, aborting; exec_report_selection')
                return
            endif
            ! read segment
            call os%new(noris)
            call bos_doc%read_segment(iseg, os)
            ! update states
            call os%set_all('state', real(states))
            call bos_doc%write_segment_inside(iseg, os)
            ! no need to update header (taken care of in binoris object)
            call bos_doc%close
        endif
    end subroutine exec_report_selection

    !> for creating a new project
    subroutine exec_new_project( self, cline )
        class(new_project_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        integer          :: iostatus
        call params%new(cline)
        if( cline%defined('projname') .and. cline%defined('dir') )then
            THROW_HARD('both projname and dir defined on command line, use either or; exec_new_project')
        else if( cline%defined('projname') )then
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

    !> for updating an existing project
    subroutine exec_update_project( self, cline )
        class(update_project_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        integer          :: iostatus
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

    !> for importing movies/integrated movies(micrographs)
    subroutine exec_import_movies( self, cline )
        class(import_movies_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        logical          :: inputted_boxtab
        integer          :: nmovf, nboxf, i
        type(ctfparams)  :: ctfvars
        character(len=:),      allocatable :: phaseplate, boxf_abspath
        character(len=LONGSTRLEN), allocatable :: boxfnames(:)
        call params%new(cline)
        ! parameter input management
        inputted_boxtab = cline%defined('boxtab')
        ! project file management
        if( .not. file_exists(trim(params%projfile)) )then
            THROW_HARD('project file: '//trim(params%projfile)//' does not exists! exec_import_movies')
        endif
        call spproj%read(params%projfile)
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
                ctfvars%ctfflag = 1
            case('no')
                ctfvars%ctfflag = 0
            case('flip')
                ctfvars%ctfflag = 2
            case DEFAULT
                THROW_HARD('ctf flag: '//trim(params%ctf)//' not supported; exec_import_movies')
        end select
        ctfvars%l_phaseplate = .false.
        if( trim(params%phaseplate) .eq. 'yes' ) ctfvars%l_phaseplate = .true.
        ! update project info
        call spproj%update_projinfo( cline )
        ! updates segment
        call spproj%add_movies(params%filetab, ctfvars)
        ! add boxtab
        if( inputted_boxtab )then
            call read_filetable(params%boxtab, boxfnames)
            nboxf = size(boxfnames)
            nmovf = nlines(params%filetab)
            if( nboxf /= nmovf )then
                write(logfhandle,*) '# boxfiles: ', nboxf
                write(logfhandle,*) '# movies  : ', nmovf
                THROW_HARD('# boxfiles .ne. # movies; exec_import_movies')
            endif
            do i=1,nmovf
                 boxf_abspath = simple_abspath(trim(boxfnames(i)), errmsg='commander_project :: exec_import_movies')
                call spproj%os_mic%set(i, 'boxfile', boxf_abspath)
            end do
        endif
        ! write project file
        call spproj%write ! full write since projinfo is updated and this is guaranteed to be the first import
        call simple_end('**** IMPORT_MOVIES NORMAL STOP ****')
    end subroutine exec_import_movies

    !> for importing particle coordinates (boxes) in EMAN1.9 format
    subroutine exec_import_boxes( self, cline )
        class(import_boxes_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        integer          :: nos_mic, nboxf, i
        character(len=:),          allocatable :: boxf_abspath
        character(len=LONGSTRLEN), allocatable :: boxfnames(:)
        call params%new(cline)
        ! project file management
        if( .not. file_exists(trim(params%projfile)) )then
            THROW_HARD('project file: '//trim(params%projfile)//' does not exist! exec_import_boxes')
        endif
        call spproj%read(params%projfile)
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
            boxf_abspath = simple_abspath(trim(boxfnames(i)),errmsg='commander_project :: exec_import_movies')
            call spproj%os_mic%set(i, 'boxfile', trim(boxf_abspath))
        end do
        ! write project file
        call spproj%write_segment_inside('mic') ! all that's needed here
        call simple_end('**** IMPORT_BOXES NORMAL STOP ****')
    end subroutine exec_import_boxes

    !> for importing extracted particles
    subroutine exec_import_particles( self, cline )
        use simple_oris,      only: oris
        use simple_nrtxtfile, only: nrtxtfile
        use simple_binoris_io ! use all in there
        class(import_particles_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        character(len=:), allocatable :: phaseplate, ctfstr
        real,             allocatable :: line(:)
        type(parameters) :: params
        type(sp_project) :: spproj
        type(oris)       :: os
        type(nrtxtfile)  :: paramfile
        logical          :: inputted_oritab, inputted_plaintexttab, inputted_deftab
        integer          :: i, ndatlines, nrecs, n_ori_inputs, lfoo(3)
        type(ctfparams)  :: ctfvars
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
        ! oris input
        if( inputted_oritab )then
            ndatlines = binread_nlines(params%oritab)
            call os%new(ndatlines)
            call binread_oritab(params%oritab, spproj, os, [1,ndatlines])
            call spproj%kill ! for safety
        endif
        if( inputted_deftab )then
            ndatlines = binread_nlines(params%deftab)
            call os%new(ndatlines)
            call binread_ctfparams_state_eo(params%deftab, spproj, os, [1,ndatlines])
            call spproj%kill ! for safety
        endif
        if( inputted_plaintexttab )then
            call paramfile%new(params%plaintexttab, 1)
            ndatlines = paramfile%get_ndatalines()
            nrecs     = paramfile%get_nrecs_per_line()
            if( nrecs < 1 .or. nrecs > 4 .or. nrecs == 2 )then
                THROW_HARD('unsupported nr of rec:s in plaintexttab; exec_extract_ptcls')
            endif
            call os%new(ndatlines)
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
                        THROW_HARD('unsupported dfunit; exec_extract_ptcls')
                end select
                select case(params%angastunit)
                    case( 'radians' )
                        if( nrecs == 3 ) line(3) = rad2deg(line(3))
                    case( 'degrees' )
                        ! nothing to do
                    case DEFAULT
                        THROW_HARD('unsupported angastunit; exec_extract_ptcls')
                end select
                select case(params%phshiftunit)
                    case( 'radians' )
                        ! nothing to do
                    case( 'degrees' )
                        if( nrecs == 4 ) line(4) = deg2rad(line(4))
                    case DEFAULT
                        THROW_HARD('unsupported phshiftunit; exec_extract_ptcls')
                end select
                call os%set(i, 'dfx', line(1))
                if( nrecs > 1 )then
                    call os%set(i, 'dfy', line(2))
                    call os%set(i, 'angast', line(3))
                endif
                if( nrecs > 3 )then
                    call os%set(i, 'phshift', line(4))
                endif
            end do
        endif
        if( cline%defined('stk') )then
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
                    THROW_HARD('unsupported ctf flag: '//trim(params%ctf)//'; exec_extract_ptcls')
            end select
            if( ctfvars%ctfflag .ne. CTFFLAG_NO )then
                ! if importing single stack of extracted particles, these are hard requirements
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
        else
            ! importing from stktab
            if( n_ori_inputs == 1 )then
                ! sampling distance
                call os%set_all2single('smpd', params%smpd)
                ! acceleration voltage
                if( cline%defined('kv') )then
                    call os%set_all2single('kv', params%kv)
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'kv') )then
                            write(logfhandle,*) 'os entry: ', i, ' lacks acceleration volatage (kv)'
                            THROW_HARD('provide kv on command line or update input document; exec_extract_ptcls')
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
                            THROW_HARD('provide cs on command line or update input document; exec_extract_ptcls')
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
                            THROW_HARD('provide fraca on command line or update input document; exec_extract_ptcls')
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

        ! PROJECT FILE MANAGEMENT
        call spproj%read(params%projfile)

        ! UPDATE FIELDS
        ! add stack if present
        if( cline%defined('stk') )then
            if( n_ori_inputs == 0 .and. trim(params%ctf) .eq. 'no' )then
                ! get number of particles from stack
                call find_ldim_nptcls(params%stk, lfoo, params%nptcls)
                call os%new(params%nptcls)
                call os%set_all2single('state', 1.0)
            endif
            call spproj%add_single_stk(params%stk, ctfvars, os)
        endif
        ! add list of stacks (stktab) if present
        if( cline%defined('stktab') ) call spproj%add_stktab(params%stktab, os)

        ! WRITE PROJECT FILE
        call spproj%write ! full write since this is guaranteed to be the first import
        call simple_end('**** IMPORT_PARTICLES NORMAL STOP ****')
    end subroutine exec_import_particles

    !> for importing class-averages projects
    subroutine exec_import_cavgs( self, cline )
        class(import_cavgs_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
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

end module simple_commander_project
