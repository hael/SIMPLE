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
    procedure :: remove
    ! writers
    procedure :: write
    ! importers
    procedure :: merge
    !other
    procedure :: kill
end type guistats

contains

    subroutine init( self, fname)
        class(guistats),            intent(inout) :: self
        character(len=*), optional, intent(in)    :: fname
        call self%stats%new(1, .false.)
        if(present(fname)) then 
            if(file_exists(fname)) call self%stats%read(fname)
        else if(file_exists(GUISTATS_FILE)) then
            call self%stats%read(GUISTATS_FILE)
        end if
    end subroutine init

    subroutine set_1( self, key, val )
        class(guistats),  intent(inout) :: self
        character(len=*), intent(in)    :: key
        real,             intent(in)    :: val
        call self%stats%set(1, key, val)
    end subroutine set_1

    subroutine set_2( self, key, val )
        class(guistats),  intent(inout) :: self
        character(len=*), intent(in)    :: key
        character(len=*), intent(in)    :: val
        call self%stats%set(1, key, val)
    end subroutine set_2

    subroutine set_3( self, key, val )
        class(guistats),  intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer,          intent(in)    :: val
        call self%stats%set(1, key, float(val))
    end subroutine set_3

    subroutine remove( self, key )
        class(guistats),  intent(inout) :: self
        character(len=*), intent(in)    :: key
        if(self%stats%isthere(1, key)) call self%stats%delete_entry(1, key)
    end subroutine remove

    subroutine merge( self, fname, delete )
        class(guistats),           intent(inout) :: self
        character(len=*),          intent(in)    :: fname
        logical,         optional, intent(in)    :: delete
        type(oris) :: instats
        type(ori)  :: instats_ori
        logical    :: l_delete
        l_delete = .false.
        if(present(delete)) l_delete = delete
        if(file_exists(fname)) then
            call instats%new(1, .false.)
            call instats%read(fname)
            if(instats%get_noris() .eq. 1) then
                call instats%get_ori(1, instats_ori)
                call self%stats%append(1, instats_ori)
                call instats_ori%kill
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
