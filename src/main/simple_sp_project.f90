module simple_sp_project
use simple_oris,    only: oris
use simple_binoris, only: binoris
use simple_fileio
implicit none

integer, parameter :: N_OS_SEG = 4

type sp_project
    ! oris representations of binary file segments
    type(oris)    :: os_stk    ! segment 1
    type(oris)    :: os_ptcl2D ! segment 2
    type(oris)    :: os_cls    ! segment 3
    type(oris)    :: os_cls3D  ! segment 4
    ! binary file-handler
    type(binoris) :: bos
contains
    procedure, private :: new_sp_oris_1
    procedure, private :: new_sp_oris_2
    generic            :: new_sp_oris => new_sp_oris_1, new_sp_oris_2
    ! I/O
    procedure          :: read
    procedure          :: write
    procedure          :: write_sp_oris
    procedure          :: read_sp_oris
    procedure          :: print_header
end type sp_project

contains

    subroutine new_sp_oris_1( self, which, os )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which
        class(oris),       intent(inout) :: os
        select case(trim(which))
            case('stk')
                self%os_stk    = os
            case('ptcl2D')
                self%os_ptcl2D = os
            case('cls')
                self%os_cls    = os
            case('cls3D')
                self%os_cls3D  = os
            case DEFAULT
                stop 'unsupported which flag; sp_project :: new_sp_oris_1'
        end select
    end subroutine new_sp_oris_1

    subroutine new_sp_oris_2( self, which, n )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which
        integer,           intent(in)    :: n
        select case(trim(which))
            case('stk')
                call self%os_stk%new_clean(n)
            case('ptcl2D')
                call self%os_ptcl2D%new_clean(n)
            case('cls')
                call self%os_cls%new_clean(n)
            case('cls3D')
                call self%os_cls3D%new_clean(n)
            case DEFAULT
                stop 'unsupported which flag; sp_project :: new_sp_oris_2'
        end select
    end subroutine new_sp_oris_2

    subroutine read( self, fname )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        integer :: isegment, n
        if( .not. file_exists(trim(fname)) )then
            write(*,*) 'fname: ', trim(fname)
            stop 'inputted file does not exist; sp_project :: new_from_file'
        endif
        if( fname2format(fname) .ne. 'O' )then
            write(*,*) 'fname: ', trim(fname)
            stop 'file format not supported; sp_project :: new_from_file'
        endif
        call self%bos%open(fname)
        do isegment=1,self%bos%get_n_segments()
            n = self%bos%get_n_records(isegment)
            select case(isegment)
                case(1)
                    call self%os_stk%new_clean(n)
                    call self%bos%read_segment(isegment, self%os_stk)
                case(2)
                    call self%os_ptcl2D%new_clean(n)
                    call self%bos%read_segment(isegment, self%os_ptcl2D)
                case(3)
                    call self%os_cls%new_clean(n)
                    call self%bos%read_segment(isegment, self%os_cls)
                case(4)
                    call self%os_cls3D%new_clean(n)
                    call self%bos%read_segment(isegment, self%os_cls3D)
            end select
        end do
        call self%bos%close
    end subroutine read

    subroutine write( self, fname )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        integer :: isegment
        if( fname2format(fname) .ne. 'O' )then
            write(*,*) 'fname: ', trim(fname)
            stop 'file format not supported; sp_project :: write'
        endif
        call self%bos%open(fname, del_if_exists=.true.)
        do isegment=1,N_OS_SEG
            select case(isegment)
                case(1)
                    if( self%os_stk%get_noris()    > 0 ) call self%bos%write_segment(isegment, self%os_stk)
                case(2)
                    if( self%os_ptcl2D%get_noris() > 0 ) call self%bos%write_segment(isegment, self%os_ptcl2D)
                case(3)
                    if( self%os_cls%get_noris()    > 0 ) call self%bos%write_segment(isegment, self%os_cls)
                case(4)
                    if( self%os_cls3D%get_noris()  > 0 ) call self%bos%write_segment(isegment, self%os_cls3D)
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
                        call self%os_stk%write(fname,    fromto)
                    case('ptcl2D')
                        call self%os_ptcl2D%write(fname, fromto)
                    case('cls')
                        call self%os_cls%write(fname,    fromto)
                    case('cls3D')
                        call self%os_cls3D%write(fname,  fromto)
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
                        call self%os_stk%read(fname,    fromto)
                    case('ptcl2D')
                        call self%os_ptcl2D%read(fname, fromto)
                    case('cls')
                        call self%os_cls%read(fname,    fromto)
                    case('cls3D')
                        call self%os_cls3D%read(fname,  fromto)
                    case DEFAULT
                        stop 'unsupported which flag; sp_project :: read_sp_oris'
                end select
            case DEFAULT
                write(*,*) 'fname: ', trim(fname)
                stop 'file format not supported; sp_project :: read_sp_oris'
        end select
    end subroutine read_sp_oris

    subroutine print_header( self )
        class(sp_project), intent(in) :: self
        call self%bos%print_header
    end subroutine print_header

end module simple_sp_project
