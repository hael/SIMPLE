!@descr: cavgs_work_project: the temporary class-average project a class-average 3D workflow refines
! Built as solve3D_cavgs builds its own: one ptcl3D row per even and per odd class average, without CTF,
! rows 1..N holding the even averages and rows N+1..2N the odd ones, every row carrying its class and the
! cls2D selection. The project registers a transient canonical sigma2 state; delete_files removes both.
module simple_cavgs_work_project
use simple_error,  only: simple_exception
use simple_oris,   only: oris
use simple_string, only: string
implicit none
private
#include "simple_local_flags.inc"

public :: cavgs_work_project

type :: cavgs_work_project
    private
    type(string) :: projfile      !< the work project file
    type(string) :: sigma2_fname  !< its transient canonical sigma2 state
    integer      :: ncavgs = 0
  contains
    procedure :: new
    procedure :: get_projfile
    procedure :: get_ncavgs
    procedure :: even_row
    procedure :: odd_row
    procedure :: cap_nparts
    procedure :: eo_angular_distance
    procedure :: eo_state_overlap
    procedure :: delete_files
    procedure :: kill
end type cavgs_work_project

contains

    !> Builds and writes the work project from the class averages spproj registers (imgkind cavg) and their
    !! _even/_odd companions; project info, compute environment and job processes are carried over
    subroutine new( self, spproj, projfile, sigma2_fname )
        use simple_fileio,     only: add2fbody, fname2ext
        use simple_ori,        only: ori
        use simple_sp_project, only: sp_project
        use simple_syslib,     only: del_file, file_exists
        use simple_type_defs,  only: CTFFLAG_NO, ctfparams
        class(cavgs_work_project), intent(inout) :: self
        class(sp_project),         intent(inout) :: spproj
        class(string),             intent(in)    :: projfile, sigma2_fname
        type(sp_project)     :: work_proj
        type(ctfparams)      :: ctfvars
        type(ori)            :: o, o_even, o_odd
        type(string)         :: stk, stk_even, stk_odd, ext
        integer, allocatable :: states(:)
        integer :: icls, even_ind, odd_ind
        real    :: smpd
        call self%kill
        call spproj%get_cavgs_stk(stk, self%ncavgs, smpd)
        if( .not. file_exists(stk) ) THROW_HARD('class-average stack does not exist: '//stk%to_char())
        ext      = string('.')//fname2ext(stk)
        stk_even = add2fbody(stk, ext, '_even')
        stk_odd  = add2fbody(stk, ext, '_odd')
        if( .not. file_exists(stk_even) ) THROW_HARD('even class-average stack does not exist: '//stk_even%to_char())
        if( .not. file_exists(stk_odd)  ) THROW_HARD('odd class-average stack does not exist: '//stk_odd%to_char())
        states = nint(spproj%os_cls2D%get_all('state'))
        if( size(states) /= self%ncavgs ) THROW_HARD('cls2D rows and registered class averages differ in number')
        if( count(states > 0) == 0 ) THROW_HARD('no selected class averages in the project')
        ctfvars%ctfflag = CTFFLAG_NO
        ctfvars%smpd    = smpd
        work_proj%projinfo = spproj%projinfo
        work_proj%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris() > 0 ) work_proj%jobproc = spproj%jobproc
        ! the rows have another identity than the input's particles: never inherit its canonical sigma2 state
        call work_proj%projinfo%delete_entry('sigma2_state')
        call work_proj%projinfo%delete_entry('projname')
        call work_proj%projinfo%delete_entry('projfile')
        call work_proj%update_projinfo(projfile)
        call del_file(sigma2_fname)
        call work_proj%set_sigma2_state_path(sigma2_fname)
        call work_proj%add_stk(stk_even, ctfvars)
        call work_proj%add_stk(stk_odd,  ctfvars)
        ! the row layout needs one even and one odd average per class
        if( work_proj%os_stk%get_int(1, 'nptcls') /= self%ncavgs ) THROW_HARD('the even stack does not hold one average per class')
        if( work_proj%os_stk%get_int(2, 'nptcls') /= self%ncavgs ) THROW_HARD('the odd stack does not hold one average per class')
        do icls = 1,self%ncavgs
            even_ind = self%even_row(icls)
            odd_ind  = self%odd_row(icls)
            call work_proj%os_ptcl3D%get_ori(icls, o)
            call o%set('class', icls)
            call o%set('state', states(icls))
            o_even = o
            call o_even%set('eo', 0)
            call o_even%set('stkind', work_proj%os_ptcl3D%get(even_ind,'stkind'))
            call work_proj%os_ptcl3D%set_ori(even_ind, o_even)
            o_odd = o
            call o_odd%set('eo', 1)
            call o_odd%set('stkind', work_proj%os_ptcl3D%get(odd_ind,'stkind'))
            call work_proj%os_ptcl3D%set_ori(odd_ind, o_odd)
            ! the selection is mirrored into ptcl2D, as everywhere else
            call work_proj%os_ptcl2D%set_state(even_ind, min(1, states(icls)))
            call work_proj%os_ptcl2D%set_state(odd_ind,  min(1, states(icls)))
        enddo
        call work_proj%write()
        self%projfile     = projfile
        self%sigma2_fname = sigma2_fname
        call work_proj%kill
        call o%kill
        call o_even%kill
        call o_odd%kill
        deallocate(states)
    end subroutine new

    function get_projfile( self ) result( fname )
        class(cavgs_work_project), intent(in) :: self
        type(string) :: fname
        fname = self%projfile
    end function get_projfile

    integer function get_ncavgs( self )
        class(cavgs_work_project), intent(in) :: self
        get_ncavgs = self%ncavgs
    end function get_ncavgs

    !> the work-project row of the even average of class icls
    integer function even_row( self, icls )
        class(cavgs_work_project), intent(in) :: self
        integer,                   intent(in) :: icls
        even_row = icls
    end function even_row

    !> the work-project row of the odd average of class icls
    integer function odd_row( self, icls )
        class(cavgs_work_project), intent(in) :: self
        integer,                   intent(in) :: icls
        odd_row = self%ncavgs + icls
    end function odd_row

    !> a distributed run gets at most one part per row
    subroutine cap_nparts( self, params, cline )
        use simple_cmdline,      only: cmdline
        use simple_defs,         only: logfhandle
        use simple_parameters,   only: parameters
        use simple_string_utils, only: int2str
        class(cavgs_work_project), intent(in)    :: self
        class(parameters),         intent(inout) :: params
        class(cmdline),            intent(inout) :: cline
        integer :: nparts_eff
        if( self%ncavgs < 1 ) THROW_HARD('cap_nparts needs a built work project')
        if( .not. cline%defined('nparts') ) return
        nparts_eff = min(params%nparts, max(1, 2 * self%ncavgs))
        if( nparts_eff < params%nparts )then
            write(logfhandle,'(A,I0,A,I0)') '>>> REDUCING NPARTS FROM ', params%nparts, &
                &' TO THE NUMBER OF EVEN/ODD CLASS AVERAGE ENTRIES: ', nparts_eff
        endif
        params%nparts = nparts_eff
        params%numlen = len(int2str(params%nparts))
        call cline%set('nparts', params%nparts)
        call cline%set('numlen', params%numlen)
    end subroutine cap_nparts

    !> the average angular distance (degrees) between the even and the odd orientation of the classes whose
    !! even row is active, over the rows of os; 0 without such classes
    real function eo_angular_distance( self, os, pgrp ) result( avg_euldist )
        use simple_linalg, only: rad2deg
        use simple_ori,    only: ori
        use simple_sym,    only: sym
        class(cavgs_work_project), intent(in) :: self
        class(oris),               intent(in) :: os
        character(len=*),          intent(in) :: pgrp
        type(sym) :: se
        type(ori) :: o_even, o_odd
        real      :: euldist
        integer   :: icls, ncls
        call se%new(pgrp)
        avg_euldist = 0.
        ncls        = 0
        do icls = 1,self%ncavgs
            call os%get_ori(self%even_row(icls), o_even)
            if( o_even%get_state() == 0 ) cycle
            ncls = ncls + 1
            call os%get_ori(self%odd_row(icls), o_odd)
            euldist = rad2deg(o_odd.euldist.o_even)
            if( se%get_nsym() > 1 )then
                call o_odd%mirror2d
                call se%rot_to_asym(o_odd)
                euldist = min(rad2deg(o_odd.euldist.o_even), euldist)
            endif
            avg_euldist = avg_euldist + euldist
        enddo
        if( ncls > 0 ) avg_euldist = avg_euldist / real(ncls)
        call o_even%kill
        call o_odd%kill
        call se%kill
    end function eo_angular_distance

    !> the percentage of classes with an active row whose even and odd rows carry the same state; 0 without
    !! such classes
    real function eo_state_overlap( self, os ) result( score )
        class(cavgs_work_project), intent(in) :: self
        class(oris),               intent(in) :: os
        integer :: icls, nsame_state, ncls, state_even, state_odd
        nsame_state = 0
        ncls        = 0
        do icls = 1,self%ncavgs
            state_even = os%get_state(self%even_row(icls))
            state_odd  = os%get_state(self%odd_row(icls))
            if( state_even == 0 .and. state_odd == 0 ) cycle
            ncls = ncls + 1
            if( state_even == state_odd ) nsame_state = nsame_state + 1
        enddo
        score = 0.
        if( ncls > 0 ) score = 100.0 * real(nsame_state) / real(ncls)
    end function eo_state_overlap

    !> removes the work project file and its sigma2 state
    subroutine delete_files( self )
        use simple_syslib, only: del_file
        class(cavgs_work_project), intent(in) :: self
        if( self%projfile%is_allocated()     ) call del_file(self%projfile)
        if( self%sigma2_fname%is_allocated() ) call del_file(self%sigma2_fname)
    end subroutine delete_files

    subroutine kill( self )
        class(cavgs_work_project), intent(inout) :: self
        call self%projfile%kill
        call self%sigma2_fname%kill
        self%ncavgs = 0
    end subroutine kill

end module simple_cavgs_work_project
