! write stats for gui 
module simple_guistats
include 'simple_lib.f08'
use simple_oris, only: oris
use simple_ori,  only: ori
implicit none

public :: guistats
private
#include "simple_local_flags.inc"

type guistats
    type(oris) :: stats
contains
    ! constructor
    procedure :: init
    ! setters 
    generic   :: set => set_1, set_2, set_3
    procedure :: set_1
    procedure :: set_2
    procedure :: set_3
    procedure :: set_now
    procedure :: remove
    ! writers
    procedure :: write
    ! importers
    procedure :: merge
    !other
    procedure :: kill
end type guistats

contains

    subroutine init( self, nlines, fname)
        class(guistats),            intent(inout) :: self
        character(len=*), optional, intent(in)    :: fname
        integer,          optional, intent(in)    :: nlines
        if(present(nlines)) then
            call self%stats%new(nlines, .false.)
        else
            call self%stats%new(1, .false.)
        end if
        if(present(fname)) then 
            if(file_exists(fname)) call self%stats%read(fname)
        else if(file_exists(GUISTATS_FILE)) then
            call self%stats%read(GUISTATS_FILE)
        end if
    end subroutine init

    subroutine set_1( self, line, key, val )
        class(guistats),  intent(inout) :: self
        character(len=*), intent(in)    :: key
        real,             intent(in)    :: val
        integer,          intent(in)    :: line
        call self%stats%set(line, key, val)
    end subroutine set_1

    subroutine set_2( self, line, key, val )
        class(guistats),  intent(inout) :: self
        character(len=*), intent(in)    :: key
        character(len=*), intent(in)    :: val
        integer,          intent(in)    :: line
        call self%stats%set(line, key, val)
    end subroutine set_2

    subroutine set_3( self, line, key, val )
        class(guistats),  intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer,          intent(in)    :: val
        integer,          intent(in)    :: line
        call self%stats%set(line, key, float(val))
    end subroutine set_3

    subroutine set_now( self, line, key )
        class(guistats),  intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer,          intent(in)    :: line
        character(8)  :: date
        character(10) :: time
        character(5)  :: zone
        character(16) :: datestr
        integer,dimension(8) :: values
        ! using keyword arguments
        call date_and_time(date,time,zone,values)
        call date_and_time(DATE=date,ZONE=zone)
        call date_and_time(TIME=time)
        call date_and_time(VALUES=values)
        write(datestr, '(I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') values(1), '/', values(2), '/', values(3), ' ', values(5), ':', values(6)
        call self%stats%set(line, key, datestr)
    end subroutine set_now

    subroutine remove( self, line, key )
        class(guistats),  intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer,          intent(in)    :: line
        if(self%stats%isthere(line, key)) call self%stats%delete_entry(line, key)
    end subroutine remove

    subroutine merge( self, fname, nlines, delete )
        class(guistats),           intent(inout) :: self
        character(len=*),          intent(in)    :: fname
        logical,         optional, intent(in)    :: delete
        integer,                   intent(in)    :: nlines
        type(oris) :: instats
        type(ori)  :: instats_ori
        logical    :: l_delete
        integer    :: line
        l_delete = .false.
        if(present(delete)) l_delete = delete
        if(file_exists(fname)) then
            call instats%new(nlines, .false.)
            call instats%read(fname)
            if(instats%get_noris() .eq. nlines) then
                do line = 1, nlines
                    call instats%get_ori(line, instats_ori)
                    call self%stats%append(line, instats_ori)
                    call instats_ori%kill
                end do
            end if
            call instats%kill
            if(l_delete) call del_file(fname)
        end if
    end subroutine merge

    subroutine write( self, fname)
        class(guistats),  intent(inout)        :: self
        character(len=*), optional, intent(in) :: fname
        if(present(fname)) then
            call self%stats%write(fname)
        else
            call self%stats%write(GUISTATS_FILE)
        end if 
    end subroutine write
    
    subroutine kill( self )
        class(guistats),  intent(inout) :: self
        call self%stats%kill
    end subroutine kill

end module simple_guistats
