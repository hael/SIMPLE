module simple_sp_project
use simple_oris,    only: oris
use simple_binoris, only: binoris
use simple_fileio
use simple_defs
implicit none

integer, parameter :: MAXN_OS_SEG = 11

type sp_project
    ! oris representations of binary file segments
    ! segments 1-10 reserved for simple program outputs, orientations and files
    type(oris)    :: os_stk    ! per-micrograph stack os, segment 1
    type(oris)    :: os_ptcl2D ! per-particle 2D os,      segment 2
    type(oris)    :: os_cls2D  ! per-cluster 2D os,       segment 3
    type(oris)    :: os_cls3D  ! per-cluster 3D os,       segment 4
    type(oris)    :: os_ptcl3D ! per-particle 3D os,      segment 5
    ! segments 11-20 reserved for project info, job management etc.
    type(oris)    :: projinfo  ! project information      segment 11
    type(oris)    :: jobproc   ! jobid + PID + etc.       segment 12
    ! binary file-handler
    type(binoris) :: bos
contains
    ! modifiers
    procedure :: new_sp_oris
    procedure :: set_sp_oris
    ! I/O
    procedure :: read
    procedure :: write
    procedure :: write_sp_oris
    procedure :: read_sp_oris
    procedure :: read_header
    procedure :: print_header
    ! destructor
    procedure :: kill
end type sp_project

contains

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
            case DEFAULT
                stop 'unsupported which flag; sp_project :: new_sp_oris_2'
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
            case DEFAULT
                stop 'unsupported which flag; sp_project :: new_sp_oris_1'
        end select
    end subroutine set_sp_oris

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
            n = self%bos%get_n_records(isegment)
            select case(isegment)
                case(STK_SEG)
                    call self%os_stk%new_clean(n)
                    call self%bos%read_segment(isegment, self%os_stk)
                case(PTCL2D_SEG)
                    call self%os_ptcl2D%new_clean(n)
                    call self%bos%read_segment(isegment, self%os_ptcl2D)
                case(CLS2D_SEG)
                    call self%os_cls2D%new_clean(n)
                    call self%bos%read_segment(isegment, self%os_cls2D)
                case(CLS3D_SEG)
                    call self%os_cls3D%new_clean(n)
                    call self%bos%read_segment(isegment, self%os_cls3D)
                case(PTCL3D_SEG)
                    call self%os_ptcl3D%new_clean(n)
                    call self%bos%read_segment(isegment, self%os_ptcl3D)
                case(PROJINFO_SEG)
                    call self%projinfo%new_clean(n)
                    call self%bos%read_segment(isegment, self%projinfo)
                case(JOBPROC_SEG)
                    call self%jobproc%new_clean(n)
                    call self%bos%read_segment(isegment, self%jobproc)
            end select
        end do
        call self%bos%close
    end subroutine read

    subroutine write( self, fname, fromto )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        integer, optional, intent(in)    :: fromto(2)
        integer :: isegment
        if( fname2format(fname) .ne. 'O' )then
            write(*,*) 'fname: ', trim(fname)
            stop 'file format not supported; sp_project :: write'
        endif
        call self%bos%open(fname, del_if_exists=.true.)
        do isegment=1,MAXN_OS_SEG
            select case(isegment)
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
                case(PROJINFO_SEG)
                    call self%bos%write_segment(isegment, self%projinfo)
                case(JOBPROC_SEG)
                    call self%bos%write_segment(isegment, self%jobproc)
            end select
        end do
        call self%bos%write_header
        call self%bos%close
    end subroutine write

    subroutine write_sp_oris( self, which, fname, fromto )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which
        character(len=*),  intent(in)    :: fname
        integer, optional, intent(in)    :: fromto(2)
        select case(fname2format(fname))
            case('O')
                write(*,*) 'fname: ', trim(fname)
                stop 'writing individual sp_oris not supported for binary *.simple files :: write_sp_oris'
            case('T')
                ! *.txt plain text ori file
                select case(trim(which))
                    case('stk')
                        if( self%os_stk%get_noris() > 0 )then
                            call self%os_stk%write(fname)
                        else
                            write(*,*) 'WARNING, no stk-type oris available to write; sp_project :: write_sp_oris'
                        endif
                    case('ptcl2D')
                        if( self%os_ptcl2D%get_noris() > 0 )then
                            call self%os_ptcl2D%write(fname, fromto)
                        else
                            write(*,*) 'WARNING, no ptcl2D-type oris available to write; sp_project :: write_sp_oris'
                        endif
                    case('cls2D')
                        if( self%os_cls2D%get_noris() > 0 )then
                            call self%os_cls2D%write(fname)
                        else
                            write(*,*) 'WARNING, no cls2D-type oris available to write; sp_project :: write_sp_oris'
                        endif
                    case('cls3D')
                        if( self%os_cls3D%get_noris() > 0 )then
                            call self%os_cls3D%write(fname,  fromto)
                        else
                            write(*,*) 'WARNING, no cls3D-type oris available to write; sp_project :: write_sp_oris'
                        endif
                    case('ptcl3D')
                        if( self%os_ptcl3D%get_noris() > 0 )then
                            call self%os_ptcl3D%write(fname, fromto)
                        else
                            write(*,*) 'WARNING, no ptcl3D-type oris available to write; sp_project :: write_sp_oris'
                        endif
                    case('projinfo')
                        if( self%projinfo%get_noris() > 0 )then
                            call self%projinfo%write(fname, fromto)
                        else
                            write(*,*) 'WARNING, no projinfo-type oris available to write; sp_project :: write_sp_oris'
                        endif
                    case('jobproc')
                        if( self%jobproc%get_noris() > 0 )then
                            call self%jobproc%write(fname)
                        else
                            write(*,*) 'WARNING, no jobproc-type oris available to write; sp_project :: write_sp_oris'
                        endif
                    case DEFAULT
                        stop 'unsupported which flag; sp_project :: write_sp_oris'
                end select
            case DEFAULT
                write(*,*) 'fname: ', trim(fname)
                stop 'file format not supported; sp_project :: write_sp_oris'
        end select
    end subroutine write_sp_oris

    subroutine read_sp_oris( self, which, fname, fromto )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which
        character(len=*),  intent(in)    :: fname
        integer, optional, intent(in)    :: fromto(2)
        select case(fname2format(fname))
            case('O')
                write(*,*) 'fname: ', trim(fname)
                stop 'reading individual sp_oris not supported for binary *.simple files :: read_sp_oris'
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
                    case('jobproc')
                        call self%jobproc%read(fname)
                    case DEFAULT
                        stop 'unsupported which flag; sp_project :: read_sp_oris'
                end select
            case DEFAULT
                write(*,*) 'fname: ', trim(fname)
                stop 'file format not supported; sp_project :: read_sp_oris'
        end select
    end subroutine read_sp_oris

    subroutine read_header( self, fname )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        if( .not. file_exists(trim(fname)) )then
            write(*,*) 'fname: ', trim(fname)
            stop 'inputted file does not exist; sp_project :: read_header'
        endif
        if( fname2format(fname) .ne. 'O' )then
            write(*,*) 'fname: ', trim(fname)
            stop 'file format not supported; sp_project :: read_header'
        endif
        call self%bos%open(fname)
        call self%bos%read_header
        call self%bos%close
    end subroutine read_header

    subroutine print_header( self )
        class(sp_project), intent(in) :: self
        call self%bos%print_header
    end subroutine print_header

    subroutine kill( self )
        class(sp_project), intent(inout) :: self
        call self%os_stk%kill
        call self%os_ptcl2D%kill
        call self%os_cls2D%kill
        call self%os_cls3D%kill
        call self%os_ptcl3D%kill
        call self%jobproc%kill
    end subroutine kill

end module simple_sp_project
