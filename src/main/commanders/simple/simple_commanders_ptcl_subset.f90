!@descr: select the best particles within balanced viewing-direction groups and export a congruent subset
module simple_commanders_ptcl_subset
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
implicit none
private
public :: commander_select_particles

type, extends(commander_base) :: commander_select_particles
  contains
    procedure :: execute => exec_select_particles
end type commander_select_particles

contains

    subroutine exec_select_particles(self, cline)
        use simple_parameters,  only: parameters
        use simple_sp_project,  only: sp_project
        use simple_ptcl_subset, only: select_ptcl_subset, write_ptcl_subset
        use simple_defs,        only: logfhandle
        use simple_string,      only: string
        use simple_error,       only: simple_exception
        use simple_syslib,      only: simple_abspath
        use simple_sym,         only: is_valid_pointgroup
        use simple_jiffys,      only: simple_end
        class(commander_select_particles), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
        class(parameters), allocatable :: params
        class(sp_project), allocatable :: project
        type(string) :: stkout, oritab, argument
        logical, allocatable :: selected(:)
        character(len=512) :: message
        integer :: status, original_count, nspace
        if( cline%defined('mkdir') )then
            if( cline%get_carg('mkdir') /= 'yes' )then
                call simple_exception('select_particles requires mkdir=yes to preserve the input project', __FILENAME__, __LINE__)
            endif
        endif
        if( .not.cline%defined('nsample') )then
            call simple_exception('select_particles requires nsample', __FILENAME__, __LINE__)
        endif
        call cline%set('mkdir','yes')
        call cline%set('oritype','ptcl3D')
        if( .not.cline%defined('nspace') ) call cline%set('nspace',64)
        nspace = cline%get_iarg('nspace')
        if( cline%get_iarg('nsample') < 1 .or. nspace < 2 .or. mod(nspace,2) /= 0 )then
            call simple_exception('select_particles requires positive nsample and a positive even nspace', __FILENAME__, __LINE__)
        endif
        if( cline%defined('state') )then
            if( cline%get_iarg('state') < 1 )then
                call simple_exception('select_particles requires a positive state', __FILENAME__, __LINE__)
            endif
        endif
        if( cline%defined('subset_score') )then
            argument = cline%get_carg('subset_score')
            if( argument /= 'corr' .and. argument /= 'corr_cart' )then
                call simple_exception('select_particles subset_score must be corr or corr_cart', __FILENAME__, __LINE__)
            endif
        endif
        if( cline%defined('pgrp') )then
            argument = cline%get_carg('pgrp')
            if( .not.is_valid_pointgroup(argument%to_char()) )then
                call simple_exception('select_particles requires a valid pgrp', __FILENAME__, __LINE__)
            endif
        endif
        if( cline%defined('eo_rebalance') )then
            argument = cline%get_carg('eo_rebalance')
            if( argument /= 'yes' .and. argument /= 'no' )then
                call simple_exception('eo_rebalance must be yes or no', __FILENAME__, __LINE__)
            endif
        endif
        allocate(params,project)
        call params%new(cline)
        call project%read(params%projfile)
        original_count = project%os_ptcl3D%get_noris()
        call select_ptcl_subset(project,params%nsample,params%nspace,trim(params%pgrp), &
            trim(params%subset_score),params%state,selected,status,message,rebalance_eo=params%eo_rebalance=='yes')
        if( status /= 0 ) call simple_exception(trim(message), __FILENAME__, __LINE__)
        if( params%eo_rebalance == 'yes' )then
            write(logfhandle,'(A)') '>>> WARNING: new direction/score-balanced particle halves split micrograph groups'
            write(logfhandle,'(A)') '>>> WARNING: original gold-standard independence is not preserved; do not reuse old half-maps/FSC'
        endif
        stkout = simple_abspath(string('selected_particles.mrc'),check_exists=.false.)
        oritab = simple_abspath(string('selected_particles.txt'),check_exists=.false.)
        call write_ptcl_subset(project,selected,stkout,oritab)
        call project%update_projinfo(params%projfile)
        call project%write(params%projfile)
        write(logfhandle,'(A,I0,A,I0)') '>>> SELECTED PARTICLES / INPUT ROWS: ',count(selected),' / ',original_count
        write(logfhandle,'(A,A)') '>>> RANKING SCORE: ',trim(params%subset_score)
        write(logfhandle,'(A,A)') '>>> SUBSET PROJECT: ',params%projfile%to_char()
        write(logfhandle,'(A,A)') '>>> SELECTED STACK: ',stkout%to_char()
        write(logfhandle,'(A,A)') '>>> ORIENTATION / CTF TABLE: ',oritab%to_char()
        call project%kill
        call simple_end('**** SIMPLE_SELECT_PARTICLES NORMAL STOP ****')
    end subroutine exec_select_particles

end module simple_commanders_ptcl_subset