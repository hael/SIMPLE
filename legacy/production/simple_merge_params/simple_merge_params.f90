!==Program simple_merge_params
!
! <merge_params/begin> is a program for merging parameters from 2D and 3D analyses. 
! <merge_params/end> 
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund X-mas 2015
!
program simple_merge_params
use simple_cmdline  ! singleton
use simple_jiffys,  only: simple_end, progress, nlines, get_fileunit, file_exists
use simple_params,  only: params
use simple_build,   only: build
use simple_oris,    only: oris
use simple_defs     ! singleton
use simple_ori,     only: ori
implicit none
type(params)            :: p
type(build)             :: b
type(oris)              :: o_clusters
type(oris), allocatable :: o_states(:)
type(ori)               :: ori3d, ori2d, ori_comp
real                    :: euls(3), corr
integer                 :: nlines_oritab, nlines_comlindoc, nlines_deftab, cnt, nexists
integer                 :: iclass, iptcl, istate, funit, nlines_state, cnt2, ntot
integer, allocatable    :: cluster2state(:), particles(:), state_translate(:,:)
logical, allocatable    :: stateset(:)
character(len=STDLEN)   :: fname
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_MERGE_PARAMS oritab=<prime2d/prime3d doc> comlindoc=<shc_clustering_nclsX.txt>'
    write(*,'(a)',advance='no') ' [doclist=<list of oritabs for the different states>] [deftab=<text file defocus values>]'
    write(*,'(a)') ' [outfile=<output parameter file{merged_params.txt}>]'
    stop
endif
call parse_cmdline
call cmdcheckvar('oritab',    1)
call cmdcheckvar('comlindoc', 2)
if( .not. defined_cmd_arg('outfile') )then
    call set_cmdline('outfile', 'merged_params.txt')
endif
call cmdcheck
p = params()                 ! parameters generated
call b%build_general_tbox(p) ! general objects built
nlines_oritab    = nlines(p%oritab)
nlines_comlindoc = nlines(p%comlindoc)
if( defined_cmd_arg('deftab') )then
    nlines_deftab = nlines(p%deftab)
    if( nlines_oritab /= nlines_deftab )then
        stop 'nr lines in oritab .ne. nr lines in deftab; must be congruent!'
    endif
endif
if( nlines_oritab  == nlines_comlindoc )then
    ! make a new oris object and read in the comlin clustering (state) info
    o_clusters = oris(nlines_comlindoc)
    call o_clusters%read(p%comlindoc)
    ! transfer state parameters to particle level
    do iptcl=1,p%nptcls
        istate = nint(o_clusters%get(iptcl, 'state'))
        call b%a%set(iptcl, 'state', real(istate))
    end do
    ! get the number of states
    p%nstates = b%a%get_nstates()
    ! get number of clusters
    p%ncls = b%a%get_ncls()
else
    p%ncls = b%a%get_ncls()
    if( p%ncls /= nlines_comlindoc )then
        stop 'nr lines in comlindoc .ne. nr clusters; must be congruent!'
    endif
    ! make a new oris object and read in the comlin clustering (state) info
    o_clusters = oris(p%ncls)
    call o_clusters%read(p%comlindoc)
    ! create mapping cluster->state
    allocate(cluster2state(p%ncls))
    do iclass=1,p%ncls
        cluster2state(iclass) = nint(o_clusters%get(iclass, 'state'))
    end do
    ! get the number of states
    p%nstates = maxval(cluster2state)
    ! transfer state parameters to particle level
    do iptcl=1,p%nptcls
        iclass = nint(b%a%get(iptcl, 'class'))
        call b%a%set(iptcl, 'state', real(cluster2state(iclass)))
    end do
endif
if( defined_cmd_arg('doclist') )then
    ! open the doclist
    funit = get_fileunit()
    open(unit=funit, status='old', file=p%doclist)
    ! read all state docs obtained with class avgs in state groups
    allocate(o_states(p%nstates), stateset(p%nptcls), state_translate(p%nstates,2))
    stateset = .false.
    nexists = 0
    do istate=1,p%nstates
        read(funit,'(a256)') fname
        if( file_exists(fname) )then
            nexists = nexists+1
            nlines_state = nlines(fname)    
            o_states(istate) = oris(nlines_state)
            call o_states(istate)%read(fname)
            state_translate(istate,1) = istate
        else
            state_translate(istate,1) = 0
        endif
    end do
    close(funit)
    if( nexists == 0 ) stop 'none of the inputted oritabs in doclist exist in the cwd'
    ntot = nexists*p%ncls
    cnt2 = 0
    do istate=1,p%nstates
        if( o_states(istate)%get_noris() > 0 )then
            cnt = 0
            do iclass=1,p%ncls
                cnt2 = cnt2+1
                call progress(cnt2, ntot)
                particles = b%a%get_cls(iclass)
                if( cluster2state(iclass) == istate )then
                    cnt = cnt+1
                    ori3d = o_states(istate)%get_ori(cnt)
                    corr  = o_states(istate)%get(cnt, 'corr')
                    do iptcl=1,size(particles)
                        ! compose ori3d and ori2d
                        ori2d    = b%a%get_ori(particles(iptcl))
                        call ori3d%compose3d2d(ori2d, ori_comp)
                        ! set parameters in b%a
                        call b%a%set_ori(particles(iptcl),ori_comp)
                        call b%a%set(particles(iptcl), 'corr',  corr)
                        stateset(particles(iptcl)) = .true.
                    end do
                endif
                deallocate(particles)
            end do
        endif
    end do
    ! wipe out stuff that was not defined in the doclists
    do iptcl=1,p%nptcls
        if( .not. stateset(iptcl) )then
            euls = 0.
            call b%a%set_euler(iptcl, euls)
            call b%a%set(iptcl, 'x',     0.)
            call b%a%set(iptcl, 'y',     0.)
            call b%a%set(iptcl, 'corr', -1.)
            call b%a%set(iptcl, 'state', 0.)
        endif
    end do
    ! re-label the states if needed
    if( nexists < p%nstates )then
        cnt = 0
        do istate=1,p%nstates
            if( state_translate(istate,1) /= 0 )then
                cnt = cnt+1
                state_translate(istate,2) = cnt
            else
                state_translate(istate,2) = 0
            endif
            if( state_translate(istate,1) /= 0 )then
                write(*,'(a,1x,i3,1x,a,1x,i3)') '>>> THE STATE THAT WAS FORMERLY:',&
                state_translate(istate,1), 'IS NOW:', state_translate(istate,2)
            else
                write(*,'(a,1x,i3,1x,a)') '>>> THE STATE THAT WAS FORMERLY:', istate, 'HAS BEEN EXCLUDED'
            endif
        end do
        do iptcl=1,p%nptcls
            istate = nint(b%a%get(iptcl, 'state'))
            call b%a%set(iptcl, 'state', real(state_translate(istate,2)))
        end do
    endif
endif
call b%a%write(p%outfile)
call simple_end('**** SIMPLE_MERGE_PARAMS NORMAL STOP ****')
end program simple_merge_params