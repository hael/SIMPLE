!@descr: select a fixed number of project micrographs into a separate non-pruned project
module simple_commanders_mic_subset
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
implicit none
private
public :: commander_select_micrographs

type, extends(commander_base) :: commander_select_micrographs
  contains
    procedure :: execute => exec_select_micrographs
end type commander_select_micrographs

contains

    subroutine exec_select_micrographs(self, cline)
        use simple_parameters, only: parameters
        use simple_sp_project, only: sp_project
        use simple_mic_subset, only: build_mic_support, rank_mic_support, select_mic_support, apply_mic_subset
        use simple_defs,       only: dp, logfhandle
        use simple_string,     only: string
        use simple_error,      only: simple_exception
        use simple_syslib,     only: file_exists, simple_abspath
        use simple_fileio,     only: fname2ext, stemname, fopen, fclose
        use simple_jiffys,     only: simple_end
        class(commander_select_micrographs), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
        class(parameters), allocatable :: params
        class(sp_project), allocatable :: project
        real(dp), allocatable :: support(:,:), before(:), after(:)
        real, allocatable :: bin_eulers(:,:)
        integer, allocatable :: populations(:), particle_mic(:)
        logical, allocatable :: eligible(:), chosen(:)
        type(string) :: output, mic_report, coverage_report, name, source_dir
        character(len=:), allocatable :: source_name
        character(len=512) :: message
        integer :: status, nodes, i, unit, bin, half, particle_target
        real(dp) :: value, upper_bound, scale
        logical :: proven, output_supplied
        output_supplied = cline%defined('outfile')
        if( cline%defined('mkdir') )then
            if( cline%get_carg('mkdir') /= 'no' )then
                call simple_exception('select_micrographs requires mkdir=no; use outfile for a separate destination', &
                                      __FILENAME__, __LINE__)
            endif
        endif
        call cline%set('mkdir','no')
        if( .not.cline%defined('nspace') ) call cline%set('nspace',64)
        allocate(params,project)
        call params%new(cline)
        if( params%prg /= 'select_micrographs' )then
            call simple_exception('select_micrographs requires explicit prg=select_micrographs', __FILENAME__, __LINE__)
        endif
        output = params%outfile
        if( .not.output_supplied )then
            source_name = trim(params%projfile%to_char())
            output = source_name(:len(source_name)-len('.simple'))//'_mic_subset.simple'
        endif
        output = simple_abspath(output,check_exists=.false.)
        if( fname2ext(output) /= 'simple' )then
            call simple_exception('outfile must have extension .simple', __FILENAME__, __LINE__)
        endif
        mic_report = output//'.micrographs.csv'
        coverage_report = output//'.coverage.csv'
        if( file_exists(output) .or. file_exists(mic_report) .or. file_exists(coverage_report) )then
            call simple_exception('subset output or report already exists; choose a new outfile', __FILENAME__, __LINE__)
        endif
        call project%read(params%projfile)
        call build_mic_support(project,trim(params%subset_mode),trim(params%subset_score),params%state, &
             params%nspace,trim(params%pgrp),support,populations,eligible,particle_mic,status,message,bin_eulers)
        if( status /= 0 ) call simple_exception(trim(message), __FILENAME__, __LINE__)
        particle_target = 0
        if( trim(params%subset_mode) == 'quality_coverage' )then
            call rank_mic_support(project,params%nmics,params%nspace,trim(params%pgrp),trim(params%subset_score), &
                params%state,populations,eligible,particle_mic,support,particle_target,status,message)
            if( status /= 0 ) call simple_exception(trim(message), __FILENAME__, __LINE__)
            write(logfhandle,'(A,I0)') '>>> DERIVED BALANCED PARTICLE TARGET: ',particle_target
        endif
        call select_mic_support(support,params%nmics,trim(params%subset_mode)=='particles', &
             params%subset_max_nodes,chosen,value,upper_bound,proven,nodes,status,candidate_mask=populations>0)
        if( status /= 0 )then
            call simple_exception('nmics exceeds eligible micrographs or selection inputs are invalid', __FILENAME__, __LINE__)
        endif
        before = sum(support,dim=2)
        after = matmul(support,merge(1._dp,0._dp,chosen))
        scale = 1._dp
        if( trim(params%subset_mode) == 'quality_coverage' )then
            scale = max(tiny(1._dp),sum(support)*real(params%nmics,dp)/ &
                real(count(populations>0),dp)/real(size(support,1),dp))
        endif
        call fopen(unit,mic_report,status='new',action='write')
        write(unit,'(A,A)') '# subset_mode=',trim(params%subset_mode)
        write(unit,'(A,A)') '# subset_score=',trim(params%subset_score)
        write(unit,'(A,I0)') '# state=',params%state
        write(unit,'(A,A)') '# pgrp=',trim(params%pgrp)
        write(unit,'(A,I0)') '# derived_particle_target=',particle_target
        write(unit,'(A,ES24.16)') '# support_scale=',scale
        write(unit,'(A,L1)') '# proven_optimal=',proven
        write(unit,'(A,I0)') '# search_nodes=',nodes
        write(unit,'(A,ES24.16)') '# objective=',value
        write(unit,'(A,ES24.16)') '# upper_bound=',upper_bound
        write(unit,'(A)') 'micrograph_index,selected,eligible_particles,selection_support,micrograph'
        do i = 1,size(chosen)
            name = ''
            if( project%os_mic%isthere(i,'intg') ) name = project%os_mic%get_str(i,'intg')
            write(unit,'(I0,A,I0,A,I0,A,ES24.16,A,A)') i,',',merge(1,0,chosen(i)),',', &
                populations(i),',',sum(support(:,i)),',',csv_field(name%to_char())
        enddo
        call fclose(unit)
        call fopen(unit,coverage_report,status='new',action='write')
        write(unit,'(A)') 'direction_bin,particle_half,e1,e2,candidate_selection_support,retained_selection_support'
        do i = 1,size(before)
            bin = (i+1)/2
            half = mod(i-1,2)
            if( trim(params%subset_mode) == 'particles' ) half = -1
            write(unit,'(I0,A,I0,A,F12.5,A,F12.5,A,ES24.16,A,ES24.16)') bin,',',half,',', &
                bin_eulers(1,bin),',',bin_eulers(2,bin),',',before(i),',',after(i)
        enddo
        call fclose(unit)
        source_dir = stemname(simple_abspath(params%projfile))
        call absolute_data_paths(project%os_stk,source_dir)
        call absolute_data_paths(project%os_mic,source_dir)
        call apply_mic_subset(project,chosen,eligible,particle_mic)
        call project%update_projinfo(output)
        call project%write(output)
        write(logfhandle,'(A,I0,A,I0)') '>>> SELECTED MICROGRAPHS / ELIGIBLE PARTICLES: ', &
            count(chosen),' / ',sum(pack(populations,chosen))
        write(logfhandle,'(A,L1,A,I0)') '>>> PROVEN OBJECTIVE OPTIMUM: ',proven,'; SEARCH NODES: ',nodes
        if( .not.proven ) write(logfhandle,'(A)') '>>> SEARCH LIMIT: subset is an incumbent, not a proven global optimum'
        write(logfhandle,'(A,A)') '>>> SUBSET PROJECT: ',output%to_char()
        call project%kill
        call simple_end('**** SIMPLE_SELECT_MICROGRAPHS NORMAL STOP ****')
    end subroutine exec_select_micrographs

    subroutine absolute_data_paths(records, source_dir)
        use simple_oris,   only: oris
        use simple_string, only: string
        use simple_syslib, only: simple_abspath
        use simple_fileio, only: filepath
        class(oris), intent(inout) :: records
        class(string), intent(in) :: source_dir
        character(len=7), parameter :: keys(4) = [character(len=7) :: 'stk','movie','intg','boxfile']
        type(string) :: path
        integer :: i, j
        do i = 1,records%get_noris()
            do j = 1,size(keys)
                if( .not.records%isthere(i,trim(keys(j))) ) cycle
                path = records%get_str(i,trim(keys(j)))
                if( path%strlen() == 0 ) cycle
                if( path%to_char([1,1]) /= '/' )then
                    path = simple_abspath(filepath(source_dir,path),check_exists=.false.)
                    call records%set(i,trim(keys(j)),path)
                endif
            enddo
        enddo
    end subroutine absolute_data_paths

    function csv_field(text) result(quoted)
        character(len=*), intent(in) :: text
        character(len=:), allocatable :: quoted
        integer :: i
        quoted = '"'
        do i = 1,len_trim(text)
            quoted = quoted//text(i:i)
            if( text(i:i) == '"' ) quoted = quoted//'"'
        enddo
        quoted = quoted//'"'
    end function csv_field

end module simple_commanders_mic_subset
