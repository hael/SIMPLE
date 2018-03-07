! concrete commander: operations on projects (spproject) and associated files
module simple_commander_project
#include "simple_lib.f08"
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_sp_project,     only: sp_project
implicit none

public :: project2txt_commander
public :: txt2project_commander
public :: print_project_info_commander
public :: manage_project_commander
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
type, extends(commander_base) :: manage_project_commander
  contains
    procedure :: execute      => exec_manage_project
end type manage_project_commander

contains

    !> convert text (.txt) oris doc to binary (.simple)
    subroutine exec_txt2project( self, cline )
        use simple_oris, only: oris
        class(txt2project_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)     :: p
        type(oris)       :: os
        type(sp_project) :: spproj
        integer          :: noris
        p     = params(cline)
        noris = nlines(p%oritab)
        call os%new_clean(noris)
        call os%read(p%oritab)
        if( file_exists(p%projfile) )then
            call spproj%read(p%projfile)
        endif
        call spproj%set_sp_oris(p%oritype, os)
        call spproj%write(p%projfile)
        call spproj%kill
        call simple_end('**** TXT2PROJECT NORMAL STOP ****')
    end subroutine exec_txt2project

    !> convert binary (.simple) oris doc to text (.txt)
    subroutine exec_project2txt( self, cline )
        class(project2txt_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)     :: p
        type(sp_project) :: spproj
        p = params(cline)
        call spproj%read(p%projfile)
        call spproj%write_segment(p%oritype, p%outfile)
        call spproj%kill
        call simple_end('**** PROJECT2TXT NORMAL STOP ****')
    end subroutine exec_project2txt

    !> convert binary (.simple) oris doc to text (.txt)
    subroutine exec_print_project_info( self, cline )
        class(print_project_info_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(params)     :: p
        type(sp_project) :: spproj
        p = params(cline)
        call spproj%read(p%projfile)
        call spproj%print_info
        call spproj%kill
        call simple_end('**** PRINT_PROJECT_INFO NORMAL STOP ****')
    end subroutine exec_print_project_info

    !> for managing projects
    subroutine exec_manage_project( self, cline )
        use simple_oris,      only: oris
        use simple_nrtxtfile, only: nrtxtfile
        use simple_binoris_io ! use all in there
        class(manage_project_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        character(len=STDLEN)         :: projfile
        character(len=:), allocatable :: phaseplate
        real,             allocatable :: line(:)
        type(sp_project) :: spproj
        type(params)     :: p
        type(oris)       :: os
        type(nrtxtfile)  :: paramfile
        logical          :: inputted_oritab, inputted_plaintexttab, inputted_deftab
        integer          :: i, ndatlines, nrecs, n_ori_inputs

        p = params(cline)

        ! PARAMETER INPUT MANAGEMENT
        ! parameter input flags
        inputted_oritab       = cline%defined('oritab')
        inputted_deftab       = cline%defined('deftab')
        inputted_plaintexttab = cline%defined('plaintexttab')
        n_ori_inputs = count([inputted_oritab,inputted_deftab,inputted_plaintexttab])
        if( n_ori_inputs > 1 )then
            write(*,*) 'ERROR, multiple parameter sources inputted, please use (oritab|deftab|plaintexttab)'
            stop 'commander_project :: exec_manage_project'
        endif
        if( cline%defined('filetab') )then
            if( n_ori_inputs > 0 )then
                write(*,*) 'Parameter input (oritab|deftab|plaintexttab) not allowed when importing movies (filetab)'
                stop 'commander_project :: exec_manage_project'
            endif
        endif
        if( inputted_oritab )then
            ndatlines = binread_nlines(p, p%oritab)
            call os%new_clean(ndatlines)
            call binread_oritab(p%oritab, spproj, os, [1,ndatlines])
            call spproj%kill ! for safety
        endif
        if( inputted_deftab )then
            ndatlines = binread_nlines(p, p%deftab)
            call os%new_clean(ndatlines)
            call binread_oritab(p%deftab, spproj, os, [1,ndatlines])
            call spproj%kill ! for safety
        endif
        if( inputted_plaintexttab )then
            call paramfile%new(p%plaintexttab, 1)
            ndatlines = paramfile%get_ndatalines()
            nrecs     = paramfile%get_nrecs_per_line()
            if( nrecs < 1 .or. nrecs > 4 .or. nrecs == 2 )then
                write(*,*) 'unsupported nr of rec:s in plaintexttab'
                stop 'simple_commander_project :: exec_manage_project'
            endif
            call os%new_clean(ndatlines)
            allocate( line(nrecs) )
            do i=1,ndatlines
                call paramfile%readNextDataLine(line)
                select case(p%dfunit)
                    case( 'A' )
                        line(1) = line(1)/1.0e4
                        if( nrecs > 1 )  line(2) = line(2)/1.0e4
                    case( 'microns' )
                        ! nothing to do
                    case DEFAULT
                        stop 'unsupported dfunit; simple_commander_project :: exec_manage_project'
                end select
                select case(p%angastunit)
                    case( 'radians' )
                        if( nrecs == 3 ) line(3) = rad2deg(line(3))
                    case( 'degrees' )
                        ! nothing to do
                    case DEFAULT
                        stop 'unsupported angastunit; simple_commander_project :: exec_manage_project'
                end select
                select case(p%phshiftunit)
                    case( 'radians' )
                        ! nothing to do
                    case( 'degrees' )
                        if( nrecs == 4 ) line(4) = deg2rad(line(4))
                    case DEFAULT
                        stop 'unsupported phshiftunit; simple_commander_project :: exec_manage_project'
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
        ! sampling distance
        if( cline%defined('smpd') )then
            call os%set_all2single('smpd', p%smpd)
        else
            do i=1,ndatlines
                if( .not. os%isthere(i, 'smpd') )then
                    write(*,*) 'os entry: ', i, ' lacks sampling distance (smpd)'
                    write(*,*) 'Please, provide smpd on command line or update input document'
                    stop 'ERROR! simple_commander_project :: exec_manage_project'
                endif
            end do
        endif
        ! acceleration voltage
        if( cline%defined('kv') )then
            call os%set_all2single('kv', p%kv)
        else
            do i=1,ndatlines
                if( .not. os%isthere(i, 'kv') )then
                    write(*,*) 'os entry: ', i, ' lacks acceleration volatage (kv)'
                    write(*,*) 'Please, provide kv on command line or update input document'
                    stop 'ERROR! simple_commander_project :: exec_manage_project'
                endif
            end do
        endif
        ! spherical aberration
        if( cline%defined('cs') )then
            call os%set_all2single('cs', p%cs)
        else
            do i=1,ndatlines
                if( .not. os%isthere(i, 'cs') )then
                    write(*,*) 'os entry: ', i, ' lacks spherical aberration constant (cs)'
                    write(*,*) 'Please, provide cs on command line or update input document'
                    stop 'ERROR! simple_commander_project :: exec_manage_project'
                endif
            end do
        endif
        ! fraction of amplitude contrast
        if( cline%defined('fraca') )then
            call os%set_all2single('fraca', p%fraca)
        else
            do i=1,ndatlines
                if( .not. os%isthere(i, 'fraca') )then
                    write(*,*) 'os entry: ', i, ' lacks fraction of amplitude contrast (fraca)'
                    write(*,*) 'Please, provide fraca on command line or update input document'
                    stop 'ERROR! simple_commander_project :: exec_manage_project'
                endif
            end do
        endif
        ! phase-plate
        if( cline%defined('phaseplate') )then
            call os%set_all2single('phaseplate', p%phaseplate)
        else
            do i=1,ndatlines
                if( .not. os%isthere(i, 'phaseplate') )then
                    call os%set(i, 'phaseplate', 'no')
                endif
            end do
        endif

        ! PROJECT FILE MANAGEMENT
        if( cline%defined('projfile') )then
            projfile = cline%get_carg('projfile')
            if( .not. file_exists(projfile) )then
                write(*,*) 'The inputted project file (projfile): ', trim(projfile)
                stop 'does not exist in cwd; commander_project :: exec_manage_project'
                call spproj%read(projfile)
            endif
        endif

        ! STACK INPUT MANAGEMENT
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            ! there needs to be associated parameters of some form
            if( n_ori_inputs < 1 )then
                write(*,*) 'ERROR, stk or stktab input requires associated parameter input (oritab|deftab|plaintexttab)'
                stop 'commander_project :: exec_manage_project'
            endif
        endif
        if( cline%defined('stk') .and. cline%defined('stktab') )then
            write(*,*) 'ERROR, stk and stktab are both defined on command line, use either or'
            stop 'commander_project :: exec_manage_project'
        endif
        if( cline%defined('filetab') )then
            if( cline%defined('stk') .or. cline%defined('stktab') )then
                write(*,*) 'ERROR, stk and stktab cannot be inputted when filetab (of movies) is inputted'
                stop 'commander_project :: exec_manage_project'
            endif
        endif

        ! UPDATE FIELDS
        ! add stack if present
        if( cline%defined('stk')     ) call spproj%add_single_stk(p%stk, os)
        ! add list of stacks (stktab) if present
        if( cline%defined('stktab')  ) call spproj%add_stktab(p%stktab, os)
        ! add list of movies (filetab) if present
        if( cline%defined('filetab') )then
            ! hard requirements
            if( .not. cline%defined('smpd')  ) stop 'smpd (sampling distance in A) input required when importing movies; commander_project :: exec_manage_project'
            if( .not. cline%defined('kv')    ) stop 'kv (acceleration voltage in kV{300}) input required when importing movies; commander_project :: exec_manage_project'
            if( .not. cline%defined('cs')    ) stop 'cs (spherical aberration constant in mm{2.7}) input required when importing movies; commander_project :: exec_manage_project'
            if( .not. cline%defined('fraca') ) stop 'fraca (fraction of amplitude contrast{0.1}) input required when importing movies; commander_project :: exec_manage_project'
            if( cline%defined('phaseplate') )then
                phaseplate = cline%get_carg('phaseplate')
            else
                allocate(phaseplate, source='no')
            endif
            call spproj%add_movies(p%filetab, p%smpd, p%kv, p%cs, p%fraca, phaseplate)
        endif
        ! update project info
        call spproj%update_projinfo( cline )
        ! update computer environment
        call spproj%update_compenv( cline )

        ! WRITE PROJECT FILE
        call spproj%write
        call simple_end('**** MANAGE_PROJECT NORMAL STOP ****')
    end subroutine exec_manage_project

end module simple_commander_project
