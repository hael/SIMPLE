!==Program simple_map2ptcls
!
! <map2ptcls/begin> is a program for mapping parameters that have been obtained using class averages to 
! the individual particle images. There are many functionalities present that will become critical in
! future releases. Right now we recommend using this program exclusively to exclude the particles
! corresponding to deselected class averages. See the workflows section of the manual for further info.
! <map2ptcls/end> 
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund X-mas 2015
!
program simple_map2ptcls
use simple_params,  only: params
use simple_build,   only: build
use simple_oris,    only: oris
use simple_ori,     only: ori
use simple_image,   only: image
use simple_cmdline, only: cmdline
use simple_corrmat
use simple_jiffys
use simple_defs
implicit none
type state_organiser
    integer, allocatable :: particles(:)
    integer   :: cls_orig = 0 
    integer   :: cls_sel  = 0
    integer   :: istate   = 0
    type(ori) :: ori3d
end type
type(params)                       :: p
type(build)                        :: b
type(cmdline)                      :: cline
type(state_organiser), allocatable :: labeler(:)
type(image), allocatable           :: imgs_sel(:), imgs_cls(:)
type(oris)                         :: o_comlindoc, o_state, a_copy, o_oritab2
type(ori)                          :: ori2d, ori_comp
real, allocatable                  :: correlations(:,:)
integer, allocatable               :: statepops(:), state_particles(:), rejected_particles(:)
integer                            :: isel, nsel, loc(1), iptcl, pind, icls
integer                            :: nlines_oritab, nlines_oritab2, nlines_comlindoc, nlines_deftab
integer                            :: cnt, istate, funit, iline, nls, lfoo(3)
logical, allocatable               :: statedoc_exists(:), selected(:)
character(len=STDLEN)              :: statedoc
real                               :: corr
logical, parameter                 :: debug=.false.
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_MAP2PTCLS stk=<particles.ext> stk2=<selected_cavgs.ext>'
    write(*,'(a)',advance='no') ' stk3=<orig_cavgs.ext> oritab=<PRIME 2D doc> [oritab2='
    write(*,'(a)',advance='no') '<prime3D shc doc>] [comlindoc=<shc_clustering_nclsX.txt>]'
    write(*,'(a)',advance='no') ' [doclist=<list of oritabs for the different states>]'
    write(*,'(a)',advance='no') ' [deftab=<text file defocus values>] [outfile=<output parameter'
    write(*,'(a)',advance='no') ' file{mapped_ptcls_params.txt}>] [mul=<shift multiplication'
    write(*,'(a)') ' factor{1}>] [nthr=<nr of OpenMP threads{1}>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',    1)
call cline%checkvar('stk2',   2)
call cline%checkvar('stk3',   3)
call cline%checkvar('oritab', 4)
if( .not. cline%defined('outfile') )then
    call cline%set('outfile', 'mapped_ptcls_params.txt')
endif
call cline%check
p = params(cline)                   ! parameters generated
call b%build_general_tbox(p, cline) ! general objects built
! find number of selected cavgs
call find_ldim_nptcls(p%stk2, lfoo, nsel)
if( debug ) print *, 'nsel: ', nsel
! find number of original cavgs
call find_ldim_nptcls(p%stk3, lfoo, p%ncls)
if( debug ) print *, 'ncls: ', p%ncls
if( p%ncls < nsel )then
    stop 'nr of original clusters cannot be less than the number of selected ones'
endif
! find number of lines in input document
nlines_oritab = nlines(p%oritab)
if( debug ) print *, 'nlines_oritab: ', nlines_oritab
if( nlines_oritab /= p%nptcls )then 
     stop 'nr lines in oritab .ne. nr images in particle stack; must be congruent!'
endif
if( cline%defined('deftab') )then
    nlines_deftab = nlines(p%deftab)
    if( nlines_oritab /= nlines_deftab )then
        stop 'nr lines in oritab .ne. nr lines in deftab; must be congruent!'
    endif
endif
if( cline%defined('doclist') )then
    if( .not. cline%defined('comlindoc') )then
        if( nlines(p%doclist) /= 1 )then
            stop 'need a comlindoc together with statelist'
        endif
    endif
endif
if( cline%defined('comlindoc') .and. cline%defined('oritab2') )then
    stop 'either comlindoc or oritab2 can be inputted, not both'
endif
if( cline%defined('doclist') .and. cline%defined('oritab2') )then
    stop 'either doclist or oritab2 can be inputted, not both'
endif
allocate(imgs_sel(nsel), imgs_cls(p%ncls))
! read images
do isel=1,nsel
    call imgs_sel(isel)%new([p%box,p%box,1], p%smpd)
    call imgs_sel(isel)%read(p%stk2, isel)
end do
do icls=1,p%ncls
    call imgs_cls(icls)%new([p%box,p%box,1], p%smpd)
    call imgs_cls(icls)%read(p%stk3, icls)
end do
write(*,'(a)') '>>> CALCULATING CORRELATIONS'
call calc_cartesian_corrmat(imgs_sel, imgs_cls, correlations)
! find selected clusters & map selected to original clusters & extract the particle indices
allocate(labeler(nsel), selected(p%ncls))
! initialise selection array
selected = .false.
write(*,'(a)') '>>> MAPPING SELECTED TO ORIGINAL CLUSTERS'
do isel=1,nsel
    loc = maxloc(correlations(isel,:))
    labeler(isel)%cls_orig  = loc(1)
    selected(loc(1)) = .true.
    if( debug ) print *, 'found orig clsind: ', labeler(isel)%cls_orig
    labeler(isel)%cls_sel   = isel
    if( debug ) print *, 'selected class index: ', labeler(isel)%cls_sel
    labeler(isel)%particles = b%a%get_cls(labeler(isel)%cls_orig)
    if( debug ) print *, 'got this number of partices: ', size(labeler(isel)%particles)
end do
! erase deselected (by setting their state to zero)
do icls=1,p%ncls
    if( selected(icls) ) cycle
    if( b%a%get_clspop(icls) > 0 )then
        rejected_particles = b%a%get_cls(icls)
        do iptcl=1,size(rejected_particles)
            call b%a%set(rejected_particles(iptcl), 'state', 0.)
        end do
        deallocate(rejected_particles)
    endif
end do
! parse state info
if( cline%defined('comlindoc') )then
    write(*,'(a)') '>>> PROCESSING COMLIN STATE ASSIGNMENT DOCUMENT'
    nlines_comlindoc = nlines(p%comlindoc)
    if( nsel /= nlines_comlindoc )then
        stop 'nr lines in comlindoc .ne. nr of selected clusters; must be congruent!'
    endif
    ! make a new oris object and read in the comlin clustering (state) info
    o_comlindoc = oris(nsel)
    call o_comlindoc%read(p%comlindoc)
    if( .not. o_comlindoc%isthere('state') )then
        write(*,*) 'no state labeling in comlindoc, perhaps you clustered with label=class'
        stop 'please, re-cluster with label=state'
    endif
    ! set the number of states
    p%nstates = o_comlindoc%get_nstates()
    ! map states to selected cavgs
    do isel=1,nsel
        labeler(isel)%istate = nint(o_comlindoc%get(isel, 'state'))
    end do
    ! extract the state populations
    allocate( statepops(p%nstates) )
    do istate=1,p%nstates
        statepops(istate) = o_comlindoc%get_statepop(istate)
    end do
else
    ! set default values for nstates, statepop and state
    p%nstates         = 1
    allocate( statepops(1) )
    statepops(1)      = nsel
    labeler(:)%istate = 1
endif
if( cline%defined('oritab2') )then
    if( .not. file_exists(p%oritab2) ) stop 'Inputted oritab2 does not exist in the cwd'
    nlines_oritab2 = nlines(p%oritab2)
    if( nlines_oritab2 /= nsel ) stop 'Nr lines in oritab2 /= nr of selected cavgs'
    o_oritab2 = oris(nsel)
    call o_oritab2%read(p%oritab2)
    ! compose orientations and set states
    do isel=1,nsel
        ! get 3d ori
        labeler(isel)%ori3d  = o_oritab2%get_ori(isel)
        labeler(isel)%istate = nint(labeler(isel)%ori3d%get('state'))
        corr                 = labeler(isel)%ori3d%get('corr')
        do iptcl=1,size(labeler(isel)%particles)
            ! get particle index 
            pind = labeler(isel)%particles(iptcl)
            ! get 2d ori
            ori2d = b%a%get_ori(pind)
            if( cline%defined('mul') )then
                call ori2d%set('x', p%mul*ori2d%get('x'))
                call ori2d%set('y', p%mul*ori2d%get('y'))
            endif
            ! transfer original parameters in b%a 
            ori_comp = b%a%get_ori(pind)
            ! compose ori3d and ori2d
            call labeler(isel)%ori3d%compose3d2d(ori2d, ori_comp)
            ! set parameters in b%a
            call b%a%set_ori(pind,ori_comp)
            call b%a%set(pind, 'corr', corr)
        end do
    end do
endif
! map states to particles
do isel=1,nsel
    do iptcl=1,size(labeler(isel)%particles)
        ! get particle index
        pind = labeler(isel)%particles(iptcl)
        call b%a%set(pind, 'state', real(labeler(isel)%istate))
    end do
end do
! parse ori info
if( cline%defined('doclist') )then
    write(*,'(a)') '>>> COMBINING 3D ORIS (CAVGS) WITH 2D ALIGNMENT (PARTICLES)'
    if( nlines(p%doclist) /= p%nstates )then
        stop 'the number of lines in doclist does not match the number of states in comlindoc'
    endif
    allocate(statedoc_exists(p%nstates))
    ! read in 3d orientations
    funit = get_fileunit()
    open(unit=funit, status='old', file=p%doclist)
    do istate=1,p%nstates
        ! read the relevant statedoc
        read(funit,'(a256)') statedoc
        statedoc_exists(istate) = file_exists(statedoc)
        if( statedoc_exists(istate) )then
            nls = nlines(statedoc)
            if( nls /= statepops(istate) )then
                write(*,*) 'the nr of lines in statedoc: ', trim(statedoc),&
                'does not match pop size: ', statepops(istate), 'in comlindoc'
                stop
            endif
            o_state = oris(nls)
            call o_state%read(statedoc)
        else
            ! make a fake o_state
            o_state = oris(statepops(istate))
            do iline=1,statepops(istate)
                call o_state%set(iline, 'state', 0.)
            end do
            statepops(istate) = 0
        endif
        cnt = 0
        do isel=1,nsel
            if( labeler(isel)%istate == istate )then
                cnt = cnt+1
                labeler(isel)%ori3d = o_state%get_ori(cnt)
            endif
        end do
    end do    
    close(funit)
    ! wipe out the states for which no docs are provided
    do isel=1,nsel
        do iptcl=1,size(labeler(isel)%particles)
            ! get particle index
            pind = labeler(isel)%particles(iptcl)
            ! get state index
            istate = nint(b%a%get(pind, 'state'))
            if( .not. statedoc_exists(istate) )then
                call b%a%set(pind, 'state', 0.)
            endif
        end do
    end do
    ! compose orientations
    do isel=1,nsel
        do iptcl=1,size(labeler(isel)%particles)
            ! get particle index
            pind = labeler(isel)%particles(iptcl)
            ! get 2d ori
            ori2d = b%a%get_ori(pind)
            if( cline%defined('mul') )then
                call ori2d%set('x', p%mul*ori2d%get('x'))
                call ori2d%set('y', p%mul*ori2d%get('y'))
            endif
            ! transfer original parameters in b%a
            ori_comp = b%a%get_ori(pind)
            ! compose ori3d and ori2d
            call labeler(isel)%ori3d%compose3d2d(ori2d, ori_comp)
            ! set parameters in b%a
            call b%a%set_ori(pind,ori_comp)
            call b%a%set(pind, 'corr',  labeler(isel)%ori3d%get('corr'))
        end do        
    end do
    ! relabel states in consequtive order
    if( any(statepops == 0) )then
        a_copy = b%a
        cnt    = 0
        do istate=1,p%nstates
            if( statepops(istate) > 0 )then
                cnt = cnt+1
                write(*,'(a,1x,i3,1x,a,1x,i3)') '>>> THE STATE THAT WAS FORMERLY:', istate, 'IS NOW:', cnt
                state_particles = nint(a_copy%get_ptcls_in_state(istate))
                do iptcl=1,size(state_particles)
                    call b%a%set(state_particles(iptcl),'state',real(cnt))
                end do
            else
                write(*,'(a,1x,i3,1x,a)') '>>> THE STATE THAT WAS FORMERLY:', istate, 'HAS BEEN EXCLUDED'
            endif
        end do
    endif
endif
call b%a%write(p%outfile)
call simple_end('**** SIMPLE_MAP2PTCLS NORMAL STOP ****')
end program simple_map2ptcls