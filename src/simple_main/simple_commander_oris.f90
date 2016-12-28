!==Class simple_commander_oris
!
! This class contains the set of concrete commanders for orientation data management in the SIMPLE library. This class provides 
! the glue between the reciver (main reciever is simple_exec program) and the abstract action, which is simply execute (defined 
! by the base class: simple_commander_base). Later we can use the composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_oris
use simple_defs            ! singleton
use simple_jiffys          ! singleton
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
implicit none

public :: cluster_oris_commander
public :: makedeftab_commander
public :: makeoris_commander
public :: map2ptcls_commander
public :: orisops_commander
public :: oristats_commander
public :: rotmats2oris_commander
private

type, extends(commander_base) :: cluster_oris_commander
 contains
   procedure :: execute      => exec_cluster_oris
end type cluster_oris_commander
type, extends(commander_base) :: makedeftab_commander
 contains
   procedure :: execute      => exec_makedeftab
end type makedeftab_commander
type, extends(commander_base) :: makeoris_commander
 contains
   procedure :: execute      => exec_makeoris
end type makeoris_commander
type, extends(commander_base) :: map2ptcls_commander
 contains
   procedure :: execute      => exec_map2ptcls
end type map2ptcls_commander
type, extends(commander_base) :: orisops_commander
 contains
   procedure :: execute      => exec_orisops
end type orisops_commander
type, extends(commander_base) :: oristats_commander
 contains
   procedure :: execute      => exec_oristats
end type oristats_commander
type, extends(commander_base) :: rotmats2oris_commander
 contains
   procedure :: execute      => exec_rotmats2oris
end type rotmats2oris_commander

contains

    subroutine exec_cluster_oris( self, cline )
        use simple_shc_cluster, only: shc_cluster
        use simple_math,        only: rad2deg
        use simple_clusterer,   only: shc_cluster_oris
        class(cluster_oris_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(build)          :: b
        type(params)         :: p
        integer              :: icls, iptcl, numlen
        real                 :: avgd, sdevd, maxd, mind
        integer, allocatable :: clsarr(:)
        p = params(cline)
        call b%build_general_tbox(p, cline)
        call shc_cluster_oris(b%a, p%ncls)
        ! calculate distance statistics
        call b%a%cluster_diststat(avgd, sdevd, maxd, mind)
        write(*,'(a,1x,f15.6)') 'AVG      GEODESIC DIST WITHIN CLUSTERS(degrees): ', rad2deg(avgd)
        write(*,'(a,1x,f15.6)') 'AVG SDEV GEODESIC DIST WITHIN CLUSTERS(degrees): ', rad2deg(sdevd)
        write(*,'(a,1x,f15.6)') 'AVG MAX  GEODESIC DIST WITHIN CLUSTERS(degrees): ', rad2deg(maxd)
        write(*,'(a,1x,f15.6)') 'AVG MIN  GEODESIC DIST WITHIN CLUSTERS(degrees): ', rad2deg(mind)
        ! generate the class documents
        numlen = len(int2str(p%ncls))
        do icls=1,p%ncls
            clsarr = b%a%get_cls_pinds(icls)
            if( allocated(clsarr) )then
                do iptcl=1,size(clsarr)
                    call b%a%write(clsarr(iptcl), 'oris_class'//int2str_pad(icls,numlen)//'.txt')
                end do
                deallocate(clsarr)
            endif
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER_ORIS NORMAL STOP ****')
    end subroutine exec_cluster_oris

    subroutine exec_makedeftab( self, cline )
        use simple_nrtxtfile, only: nrtxtfile
        use simple_math,      only: rad2deg
        class(makedeftab_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(build)       :: b
        type(params)      :: p
        integer           :: nptcls, iptcl, ndatlines, nrecs
        type(nrtxtfile)   :: ctfparamfile
        real, allocatable :: line(:)
        p = params(cline)
        call b%build_general_tbox(p, cline)
        if( cline%defined('oritab') .or. cline%defined('deftab')  )then
            if( cline%defined('oritab') )then
                nptcls = nlines(p%oritab)
                call b%a%new(nptcls)
                call b%a%read(p%oritab)
            endif
            if( cline%defined('deftab') )then
                nptcls = nlines(p%deftab)
                call b%a%new(nptcls)
                call b%a%read(p%deftab)
            endif
            if( b%a%isthere('dfx') .and. b%a%isthere('dfy') .and. b%a%isthere('angast') )then
                ! all ok, astigmatic CTF
            else if( b%a%isthere('dfx') )then
                ! all ok, non-astigmatic CTF
            else
                write(*,*) 'defocus params (dfx, dfy, anagst) need to be in inputted via oritab/deftab'
                stop 'simple_commander_oris :: exec_makedeftab'
            endif
            do iptcl=1,nptcls
                call b%a%set(iptcl, 'kv',    p%kv   )
                call b%a%set(iptcl, 'cs',    p%cs   )
                call b%a%set(iptcl, 'fraca', p%fraca)
            end do
        else if( cline%defined('plaintexttab') )then
            call ctfparamfile%new(p%plaintexttab, 1)
            ndatlines = ctfparamfile%get_ndatalines()
            nrecs     = ctfparamfile%get_nrecs_per_line()
            if( nrecs < 1 .or. nrecs > 3 .or. nrecs == 2 )then
                write(*,*) 'unsupported nr of rec:s in plaintexttab'
                stop 'simple_commander_oris :: exec_makedeftab'
            endif
            call b%a%new(ndatlines)
            allocate( line(nrecs) )
            do iptcl=1,ndatlines
                call ctfparamfile%readNextDataLine(line)
                select case(p%dfunit)
                    case( 'A' )
                        line(1) = line(1)/1.0e4
                        if( nrecs > 1 )  line(2) = line(2)/1.0e4
                    case( 'microns' )
                        ! nothing to do
                    case DEFAULT 
                        stop 'unsupported dfunit; simple_commander_oris :: exec_makedeftab'
                end select
                select case(p%angastunit)
                    case( 'radians' )
                        if( nrecs == 3 ) line(3) = rad2deg(line(3))
                    case( 'degrees' )
                        ! nothing to do
                    case DEFAULT 
                        stop 'unsupported angastunit; simple_commander_oris :: exec_makedeftab'
                end select
                call b%a%set(iptcl, 'kv',    p%kv   )
                call b%a%set(iptcl, 'cs',    p%cs   )
                call b%a%set(iptcl, 'fraca', p%fraca)
                call b%a%set(iptcl, 'dfx',   line(1))
                if( nrecs > 1 )then
                    call b%a%set(iptcl, 'dfy',    line(2))
                    call b%a%set(iptcl, 'angast', line(3))
                endif
            end do
        else
            write(*,*) 'Nothing to do!'
        endif
        call b%a%write(p%outfile)
        ! end gracefully
        call simple_end('**** SIMPLE_MAKEDEFTAB NORMAL STOP ****')
    end subroutine exec_makedeftab

    subroutine exec_makeoris( self, cline )
        use simple_ori,  only: ori
        use simple_oris, only: oris
        use simple_math, only: normvec
        use simple_rnd,  only: irnd_uni, ran3
        class(makeoris_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(build)  :: b
        type(ori)    :: orientation
        type(oris)   :: o, o_even
        type(params) :: p
        real         :: e3, x, y
        integer      :: i, j, cnt, ispace, irot, class
        p = params(cline)
        call b%build_general_tbox(p, cline)
        if( cline%defined('ncls') )then
            if( cline%defined('angerr') )then
                o_even = oris(p%ncls)
                call o_even%spiral(p%nsym, p%eullims)
                call b%a%new(p%nptcls)
                do i=1,p%nptcls
                    class = irnd_uni(p%ncls)
                    orientation = o_even%get_ori(class)
                    call b%a%set_ori(i, orientation)
                    e3 = ran3()*2.*p%angerr-p%angerr
                    x  = ran3()*2.0*p%trs-p%trs
                    y  = ran3()*2.0*p%trs-p%trs
                    call b%a%set(i, 'x', x)
                    call b%a%set(i, 'y', y)
                    call b%a%e3set(i, e3)
                    call b%a%set(i, 'class', real(class))
                end do
            else
                o = oris(p%ncls)
                call o%spiral(p%nsym, p%eullims)
                call b%a%new(p%ncls*p%minp)
                cnt = 0
                do i=1,p%ncls
                    orientation = o%get_ori(i)
                    do j=1,p%minp
                        cnt = cnt+1
                        call b%a%set_ori(cnt, orientation)
                    end do
                end do
                if( p%zero .ne. 'yes' ) call b%a%rnd_inpls(p%trs) 
            endif
        else if( cline%defined('ndiscrete') )then
            if( p%ndiscrete > 0 )then
                call b%a%rnd_oris_discrete(p%ndiscrete, p%nsym, p%eullims)
            endif
            call b%a%rnd_inpls(p%trs)
        else if( p%even .eq. 'yes' )then
            call b%a%spiral(p%nsym, p%eullims)
            call b%a%rnd_inpls(p%trs)
        else if( p%diverse .eq. 'yes' )then
            call b%a%gen_diverse
        else if( cline%defined('nspace') )then
            ! create the projection directions
            call o%new(p%nspace)
            call o%spiral
            ! count the number of orientations
            cnt = 0
            do ispace=1,p%nspace
                do irot=0,359,p%iares
                    cnt = cnt+1
                end do
            end do
            ! fill up
            call b%a%new(cnt)
            cnt = 0
            do ispace=1,p%nspace
                orientation = o%get_ori(ispace)
                 do irot=0,359,p%iares
                    cnt = cnt+1
                    call orientation%e3set(real(irot))
                    call b%a%set_ori(cnt, orientation)
                end do
            end do
        else
            call b%a%rnd_oris(p%trs) 
        endif
        if( p%nstates > 1 ) call b%a%rnd_states(p%nstates)
        if( cline%defined('astigerr') )then
            if( p%ctf .eq. 'yes' ) call b%a%rnd_ctf(p%kv, p%cs, p%fraca, p%defocus, p%dferr, p%astigerr)
        else
            if( p%ctf .eq. 'yes' ) call b%a%rnd_ctf(p%kv, p%cs, p%fraca, p%defocus, p%dferr)
        endif
        call b%a%write(p%outfile)
        ! end gracefully
        call simple_end('**** SIMPLE_MAKEORIS NORMAL STOP ****')
    end subroutine exec_makeoris
    
    subroutine exec_map2ptcls( self, cline )
        use simple_oris,    only: oris
        use simple_ori,     only: ori
        use simple_image,   only: image
        use simple_corrmat  ! singleton
        class(map2ptcls_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type state_organiser
            integer, allocatable :: particles(:)
            integer              :: cls_orig = 0 
            integer              :: cls_sel  = 0
            integer              :: istate   = 0
            type(ori)            :: ori3d
        end type state_organiser
        type(params)                       :: p
        type(build)                        :: b
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
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        ! find number of selected cavgs
        call find_ldim_nptcls(p%stk2, lfoo, nsel)
        if( debug ) print *, 'nsel: ', nsel
        ! find number of original cavgs
        call find_ldim_nptcls(p%stk3, lfoo, p%ncls)
        if( debug ) print *, 'ncls: ', p%ncls
        if( p%ncls < nsel ) stop 'nr of original clusters cannot be less than the number of selected ones'
        ! find number of lines in input document
        nlines_oritab = nlines(p%oritab)
        if( debug ) print *, 'nlines_oritab: ', nlines_oritab
        if( nlines_oritab /= p%nptcls ) stop 'nr lines in oritab .ne. nr images in particle stack; must be congruent!'
        if( cline%defined('deftab') )then
            nlines_deftab = nlines(p%deftab)
            if( nlines_oritab /= nlines_deftab ) stop 'nr lines in oritab .ne. nr lines in deftab; must be congruent!'
        endif
        if( cline%defined('doclist') )then
            if( .not. cline%defined('comlindoc') )then
                if( nlines(p%doclist) /= 1 ) stop 'need a comlindoc together with statelist'        
            endif
        endif
        if( cline%defined('comlindoc') .and. cline%defined('oritab2') ) stop 'either comlindoc or oritab2 can be inputted, not both'
        if( cline%defined('doclist')   .and. cline%defined('oritab2') ) stop 'either doclist or oritab2 can be inputted, not both'
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
            loc                     = maxloc(correlations(isel,:))
            labeler(isel)%cls_orig  = loc(1)
            selected(loc(1))        = .true.
            if( debug ) print *, 'found orig clsind: ', labeler(isel)%cls_orig
            labeler(isel)%cls_sel   = isel
            if( debug ) print *, 'selected class index: ', labeler(isel)%cls_sel
            labeler(isel)%particles = b%a%get_cls_pinds(labeler(isel)%cls_orig)
            if( debug ) print *, 'got this number of partices: ', size(labeler(isel)%particles)
        end do
        ! erase deselected (by setting their state to zero)
        do icls=1,p%ncls
            if( selected(icls) ) cycle
            if( b%a%get_clspop(icls) > 0 )then
                rejected_particles = b%a%get_cls_pinds(icls)
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
            if( nsel /= nlines_comlindoc ) stop 'nr lines in comlindoc .ne. nr of selected clusters; must be congruent!'
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
    end subroutine exec_map2ptcls
   
    subroutine exec_orisops(self,cline)
        use simple_ori,  only: ori
        use simple_math, only: normvec
        class(orisops_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(build)  :: b
        type(ori)    :: orientation
        type(params) :: p
        real         :: normal(3)
        integer      :: s, i
        p = params(cline)
        call b%build_general_tbox(p, cline)
        if( p%errify .eq. 'yes' )then   ! introduce error in input orientations
            if( cline%defined('angerr') .or. cline%defined('sherr') ) call b%a%introd_alig_err(p%angerr, p%sherr)
            if( p%ctf .eq. 'yes' ) call b%a%introd_ctf_err(p%dferr)
        endif
        if( p%mirr .eq. '2d' ) call b%a%mirror2d ! mirror input Eulers
        if( p%mirr .eq. '3d' ) call b%a%mirror3d ! mirror input Eulers
        if( cline%defined('e1') )then ! rotate input Eulers
            call orientation%new
            call orientation%e1set(p%e1)
            call orientation%e2set(p%e2)
            call orientation%e3set(p%e3) 
            if( cline%defined('state') )then
                do i=1,b%a%get_noris()
                    s = nint(b%a%get(i, 'state'))
                    if( s == p%state )then
                        call b%a%rot(i,orientation)
                    endif
                end do
            else
                call b%a%rot(orientation)
            endif
        endif
        if( cline%defined('mul') )then
            call b%a%mul_shifts(p%mul)
        endif
        if( p%zero  .eq. 'yes' ) call b%a%zero_shifts
        if( p%plot  .eq. 'yes' )then ! plot polar vectors                          
            do i=1,b%a%get_noris()
                normal = b%a%get_normal(i)
                write(*,'(1x,f7.2,3x,f7.2)') normal(1), normal(2)
            end do
        endif
        if( p%discrete .eq. 'yes' )then
            if( cline%defined('ndiscrete') )then
                call b%a%discretize(p%ndiscrete)
            else
                stop 'need ndiscrete to be defined!'
            endif
        endif
        if( cline%defined('xsh') )     call b%a%map3dshift22d([p%xsh,p%ysh,p%zsh])
        if( cline%defined('nstates') ) call b%a%rnd_states(p%nstates)
        call b%a%write(p%outfile)
        call simple_end('**** SIMPLE_ORISOPS NORMAL STOP ****')
    end subroutine exec_orisops
    
    subroutine exec_oristats(self,cline)
        use simple_oris, only: oris
        class(oristats_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(build)  :: b
        type(oris)   :: o
        type(params) :: p
        real         :: mind, maxd, avgd, sdevd, sumd, vard
        real         :: mind2, maxd2, avgd2, sdevd2, vard2, homo_cnt, homo_avg
        logical      :: err
        p = params(cline)
        call b%build_general_tbox(p, cline)
        if( cline%defined('oritab2') )then
            ! Comparison
            if( .not. cline%defined('oritab') ) stop 'need oritab for comparison'
            if( nlines(p%oritab) .ne. nlines(p%oritab2) )then
                stop 'inconsistent number of lines in the two oritabs!'
            endif
            o = oris(p%nptcls)
            call o%read(p%oritab2)
            call b%a%diststat(o, sumd, avgd, sdevd, mind, maxd)
            write(*,'(a,1x,f15.6)') 'SUM OF ANGULAR DISTANCE BETWEEN ORIENTATIONS  :', sumd
            write(*,'(a,1x,f15.6)') 'AVERAGE ANGULAR DISTANCE BETWEEN ORIENTATIONS :', avgd
            write(*,'(a,1x,f15.6)') 'STANDARD DEVIATION OF ANGULAR DISTANCES       :', sdevd
            write(*,'(a,1x,f15.6)') 'MINIMUM ANGULAR DISTANCE                      :', mind
            write(*,'(a,1x,f15.6)') 'MAXIMUM ANGULAR DISTANCE                      :', maxd
        else if( cline%defined('oritab') )then
            ! General info
            if( cline%defined('hist') )then
                call b%a%histogram(p%hist)
                goto 999
            endif
            if( p%ctfstats .eq. 'yes' )then
                call b%a%stats('ctfres', avgd, sdevd, vard, err )
                call b%a%minmax('ctfres', mind, maxd)
                write(*,'(a,1x,f8.2)') 'AVERAGE CTF RESOLUTION               :', avgd
                write(*,'(a,1x,f8.2)') 'STANDARD DEVIATION OF CTF RESOLUTION :', sdevd
                write(*,'(a,1x,f8.2)') 'MINIMUM CTF RESOLUTION (BEST)        :', mind
                write(*,'(a,1x,f8.2)') 'MAXIMUM CTF RESOLUTION (WORST)       :', maxd
                call b%a%stats('dfx', avgd, sdevd, vard, err )
                call b%a%minmax('dfx', mind, maxd)
                call b%a%stats('dfy', avgd2, sdevd2, vard2, err )
                call b%a%minmax('dfy', mind2, maxd2)
                write(*,'(a,1x,f8.2)') 'AVERAGE DF                           :', (avgd+avgd2)/2.
                write(*,'(a,1x,f8.2)') 'STANDARD DEVIATION OF DF             :', (sdevd+sdevd2)/2.
                write(*,'(a,1x,f8.2)') 'MINIMUM DF                           :', (mind+mind2)/2.
                write(*,'(a,1x,f8.2)') 'MAXIMUM DF                           :', (maxd+maxd2)/2.
                goto 999
            endif
            if( p%trsstats .eq. 'yes' )then
                call b%a%stats('x', avgd, sdevd, vard, err )
                call b%a%minmax('x', mind, maxd)
                call b%a%stats('y', avgd2, sdevd2, vard2, err )
                call b%a%minmax('y', mind2, maxd2)
                write(*,'(a,1x,f8.2)') 'AVERAGE TRS               :', (avgd+avgd2)/2.
                write(*,'(a,1x,f8.2)') 'STANDARD DEVIATION OF TRS :', (sdevd+sdevd2)/2.
                write(*,'(a,1x,f8.2)') 'MINIMUM TRS               :', (mind+mind2)/2.
                write(*,'(a,1x,f8.2)') 'MAXIMUM TRS               :', (maxd+maxd2)/2.
                goto 999
            endif
            ! Class and states
            if( p%clustvalid .eq. 'yes' )then
                if( cline%defined('ncls') )then
                    write(*,'(a,3x,f5.1)') '>>> COHESION: ',   b%a%cohesion_norm('class',p%ncls)*100.
                    write(*,'(a,1x,f5.1)') '>>> SEPARATION: ', b%a%separation_norm('class',p%ncls)*100.
                else if( cline%defined('nstates') )then
                    write(*,'(a,3x,f5.1)') '>>> COHESION: ',   b%a%cohesion_norm('state',p%nstates)*100.
                    write(*,'(a,1x,f5.1)') '>>> SEPARATION: ', b%a%separation_norm('state',p%nstates)*100.
                else
                    stop 'need ncls/nstates as input for clustvalid'
                endif
            else if( p%clustvalid .eq. 'homo' )then
                if( cline%defined('ncls') )then
                    call b%a%homogeneity('class', p%minp, p%thres, homo_cnt, homo_avg)
                    write(*,'(a,1x,f5.1)') '>>> THIS % OF CLUSTERS CONSIDERED HOMOGENEOUS: ', homo_cnt*100.
                    write(*,'(a,1x,f5.1)') '>>> AVERAGE HOMOGENEITY:                       ', homo_avg*100.
                else if( cline%defined('nstates') )then
                    call b%a%homogeneity('state', p%minp, p%thres, homo_cnt, homo_avg)
                    write(*,'(a,1x,f5.1)') '>>> THIS % OF CLUSTERS CONSIDERED HOMOGENEOUS: ', homo_cnt*100.
                    write(*,'(a,13x,f5.1)') '>>> AVERAGE HOMOGENEITY:                      ', homo_avg*100.
                else
                    stop 'need ncls/nstates as input for clustvalid'
                endif
            endif
        endif
        call b%a%write(p%outfile)
        999 call simple_end('**** SIMPLE_ORISTATS NORMAL STOP ****')
    end subroutine exec_oristats

    subroutine exec_rotmats2oris( self, cline )
        use simple_oris,      only: oris
        use simple_ori,       only: ori
        use simple_nrtxtfile, only: nrtxtfile
        class(rotmats2oris_commander),  intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(nrtxtfile) :: rotmats
        type(params)    :: p
        type(ori)       :: o
        type(oris)      :: os_out
        integer         :: nrecs_per_line, iline, ndatlines
        real            :: rline(9), rmat(3,3)
        p = params(cline)
        if( cline%defined('infile') )then
            call rotmats%new(p%infile, 1)
            ndatlines = rotmats%get_ndatalines()
        else
            stop 'Need infile defined on command line: text file with 9 &
            &records per line defining a rotation matrix (11) (12) (13) (21) etc.'
        endif
        nrecs_per_line = rotmats%get_nrecs_per_line()
        if( nrecs_per_line /= 9 ) stop 'need 9 records (real nrs) per&
        &line of file (infile) describing rotation matrices'
        call os_out%new(ndatlines)
        do iline=1,ndatlines
            call rotmats%readNextDataLine(rline)
            ! print *, rline(1), rline(2), rline(3), rline(4), rline(5), rline(6), rline(7), rline(8), rline(9)
            rmat(1,1) = rline(1)
            rmat(1,2) = rline(2)
            rmat(1,3) = rline(3)
            rmat(2,1) = rline(4)
            rmat(2,2) = rline(5)
            rmat(2,3) = rline(6)
            rmat(3,1) = rline(7)
            rmat(3,2) = rline(8)
            rmat(3,3) = rline(9)
            rmat = transpose(rmat)
            call o%ori_from_rotmat(rmat)
            call os_out%set_ori(iline,o)
        end do
        call os_out%swape1e3
        call os_out%write(p%outfile)
        call rotmats%kill
        call simple_end('**** ROTMATS2ORIS NORMAL STOP ****')
    end subroutine exec_rotmats2oris

end module simple_commander_oris
