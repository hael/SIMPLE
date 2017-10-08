! concrete commander: operations on orientations
module simple_commander_oris
#include "simple_lib.f08"

use simple_ori,            only: ori
use simple_oris,           only: oris
use simple_binoris_io,     only: binwrite_oritab, binread_oritab,binread_nlines
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_nrtxtfile,      only: nrtxtfile
implicit none

public :: cluster_oris_commander
public :: makedeftab_commander
public :: makeoris_commander
public :: map2ptcls_commander
public :: orisops_commander
public :: oristats_commander
public :: rotmats2oris_commander
public :: bin2txt_commander
public :: txt2bin_commander
public :: vizoris_commander

private
#include "simple_local_flags.inc"

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
type, extends(commander_base) :: bin2txt_commander
 contains
   procedure :: execute      => exec_bin2txt
end type bin2txt_commander
type, extends(commander_base) :: txt2bin_commander
 contains
   procedure :: execute      => exec_txt2bin
end type txt2bin_commander
type, extends(commander_base) :: vizoris_commander
 contains
   procedure :: execute      => exec_vizoris
end type vizoris_commander

contains

    !> cluster_oris is a program for clustering orientations based on geodesic distance
    subroutine exec_cluster_oris( self, cline )
        use simple_shc_cluster, only: shc_cluster
        use simple_clusterer,   only: shc_cluster_oris
        use simple_oris,        only: oris
        class(cluster_oris_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(build)          :: b
        type(params)         :: p
        type(oris)           :: os_class
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
            clsarr = b%a%get_pinds(icls, 'class')                               !! realloc warning
            if( allocated(clsarr) )then
                call os_class%new(size(clsarr))
                do iptcl=1,size(clsarr)
                    call os_class%set_ori(iptcl, b%a%get_ori(clsarr(iptcl)))
                end do 
                call binwrite_oritab('oris_class'//int2str_pad(icls,numlen)//METADATEXT,&
                    &os_class, [1,size(clsarr)])
                deallocate(clsarr)
                call os_class%kill
            endif
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER_ORIS NORMAL STOP ****')
    end subroutine exec_cluster_oris

    !>  for creating a SIMPLE conformant file of CTF parameter values (deftab)
    subroutine exec_makedeftab( self, cline )
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
                nptcls = binread_nlines(p%oritab)
                call b%a%new(nptcls)
                call binread_oritab(p%oritab, b%a, [1,nptcls])
            endif
            if( cline%defined('deftab') )then
                nptcls = binread_nlines(p%deftab)
                call b%a%new(nptcls)
                call binread_oritab(p%deftab, b%a, [1,nptcls])
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
                call b%a%set(iptcl, 'smpd',  p%smpd )
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
                call b%a%set(iptcl, 'smpd',  p%smpd )
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
        call binwrite_oritab(p%outfile, b%a, [1,b%a%get_noris()])
        ! end gracefully
        call simple_end('**** SIMPLE_MAKEDEFTAB NORMAL STOP ****')
    end subroutine exec_makedeftab

    !> for making SIMPLE orientation/parameter files
    subroutine exec_makeoris( self, cline )
        use simple_ori,           only: ori
        use simple_oris,          only: oris
        use simple_combinatorics, only: shc_aggregation
        class(makeoris_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        character(len=STDLEN), allocatable :: oritabs(:)
        integer, allocatable :: labels(:,:), consensus(:)
        real,    allocatable :: corrs(:)
        type(build)  :: b
        type(ori)    :: orientation
        type(oris)   :: o, o_even
        type(params) :: p
        real         :: e3, x, y!, score
        integer      :: i, j, cnt, ispace, irot, class, loc(1)
        integer      :: ioritab, noritabs, nl, nl1
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
        else if( cline%defined('doclist') )then
            call read_filetable(p%doclist, oritabs)
            noritabs = size(oritabs)
            do ioritab=1,noritabs
                nl = binread_nlines(oritabs(ioritab))
                if( ioritab == 1 )then
                    nl1 = nl
                else
                    if( nl /= nl1 ) stop 'nonconfoming nr of oris in oritabs;&
                    &simple_commander_oris :: makeoris'
                endif
            end do
            allocate(labels(noritabs,nl), consensus(nl), corrs(noritabs))
            call o%new(nl)
            do ioritab=1,noritabs
                call binread_oritab(oritabs(ioritab), o, [1,nl])
                labels(ioritab,:) = nint(o%get_all('state'))
                corrs(ioritab)    = sum(o%get_all('corrs'), mask=(labels(ioritab,:)>0))&
                &/real(count(labels(ioritab,:)>0))
            end do
            loc = maxloc(corrs)
            call shc_aggregation(noritabs, nl, labels, consensus, loc(1))
            do i=1,nl
                call o%set(i,'state', real(consensus(i)))
            end do
            call binwrite_oritab('aggregated_oris'//METADATEXT, o, [1,nl])
            return
        else
            call b%a%rnd_oris(p%trs)
            if( p%doprint .eq. 'yes' )then
                call b%a%print_matrices
            endif
        endif
        if( p%nstates > 1 ) call b%a%rnd_states(p%nstates)
        if( cline%defined('astigerr') )then
            if( p%ctf .eq. 'yes' ) call b%a%rnd_ctf(p%kv, p%cs, p%fraca, p%defocus, p%dferr, p%astigerr)
        else
            if( p%ctf .eq. 'yes' ) call b%a%rnd_ctf(p%kv, p%cs, p%fraca, p%defocus, p%dferr)
        endif
        call binwrite_oritab(p%outfile, b%a, [1,b%a%get_noris()])
        ! end gracefully
        call simple_end('**** SIMPLE_MAKEORIS NORMAL STOP ****')
    end subroutine exec_makeoris

    !> for mapping parameters that have been obtained using class averages to the individual particle images
    subroutine exec_map2ptcls( self, cline )
        use simple_oris,    only: oris
        use simple_ori,     only: ori
        use simple_image,   only: image
        use simple_imghead,  only: find_ldim_nptcls
        use simple_corrmat   ! use all in there
        class(map2ptcls_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type state_organiser !> map2ptcls state struct
            integer, allocatable :: particles(:)
            type(ori)            :: ori3d
        end type state_organiser
        type(state_organiser), allocatable :: labeler(:)
        type(image),           allocatable :: imgs_sel(:), imgs_cls(:)
        real,                  allocatable :: correlations(:,:)
        integer,               allocatable :: rejected_particles(:)
        logical,               allocatable :: selected(:)
        integer      :: isel, nsel, loc(1), iptcl, pind, icls
        integer      :: nlines_oritab, nlines_oritab3D, nlines_deftab
        integer      :: lfoo(3)
        real         :: corr, rproj, rstate
        type(params) :: p
        type(build)  :: b
        type(oris)   :: o_oritab3D
        type(ori)    :: ori2d, ori_comp, o
        if( cline%defined('doclist')   ) stop 'doclist execution route no longer supported'
        if( cline%defined('comlindoc') ) stop 'comlindoc execution route no longer supported'
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        ! find number of selected cavgs
        call find_ldim_nptcls(p%stk2, lfoo, nsel)
        ! find number of original cavgs
        call find_ldim_nptcls(p%stk3, lfoo, p%ncls)
        if( p%ncls < nsel ) stop 'nr of original clusters cannot be less than the number of selected ones'
        ! find number of lines in input document
        nlines_oritab = binread_nlines(p%oritab)
        if( nlines_oritab /= p%nptcls ) stop 'nr lines in oritab .ne. nr images in particle stack; must be congruent!'
        if( cline%defined('deftab') )then
            nlines_deftab = binread_nlines(p%deftab)
            if( nlines_oritab /= nlines_deftab ) stop 'nr lines in oritab .ne. nr lines in deftab; must be congruent!'
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
            loc                     = maxloc(correlations(isel,:))
            selected(loc(1))        = .true.
            labeler(isel)%particles = b%a%get_pinds(loc(1), 'class')
        end do
        ! erase deselected (by setting their state to zero)
        do icls=1,p%ncls
            if( selected(icls) ) cycle
            if( b%a%get_pop(icls, 'class') > 0 )then
                rejected_particles = b%a%get_pinds(icls, 'class')
                do iptcl=1,size(rejected_particles)
                    call b%a%set(rejected_particles(iptcl), 'state', 0.)
                end do
                deallocate(rejected_particles)
            endif
        end do
        if( cline%defined('oritab3D') )then
            if( .not. file_exists(p%oritab3D) ) stop 'Inputted oritab3D does not exist in the cwd'
            nlines_oritab3D = binread_nlines(p%oritab3D)
            if( nlines_oritab3D /= nsel ) stop '# lines in oritab3D /= nr of selected cavgs'
            o_oritab3D = oris(nsel)
            call binread_oritab(p%oritab3D, o_oritab3D, [1,nsel])
            ! compose orientations and set states
            do isel=1,nsel
                ! get 3d ori info
                o      = o_oritab3D%get_ori(isel)
                rproj  = o%get('proj')
                rstate = o%get('state')
                corr   = o%get('corr')
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
                    call o%compose3d2d(ori2d, ori_comp)
                    ! set parameters in b%a
                    call b%a%set_ori(pind, ori_comp)
                    call b%a%set(pind, 'corr',  corr)
                    call b%a%set(pind, 'proj',  rproj)
                    call b%a%set(pind, 'state', rstate)
                end do
            end do
        endif
        call binwrite_oritab(p%outfile, b%a, [1,p%nptcls])
        call simple_end('**** SIMPLE_MAP2PTCLS NORMAL STOP ****')
    end subroutine exec_map2ptcls

    subroutine exec_orisops(self,cline)
        use simple_ori,  only: ori
        use simple_math, only: normvec
        use simple_math, only: hpsort
        class(orisops_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(build)  :: b
        type(ori)    :: orientation
        type(params) :: p
        integer      :: s, i
        p = params(cline)
        call b%build_general_tbox(p, cline)
        if( p%errify .eq. 'yes' )then   ! introduce error in input orientations
            if( cline%defined('angerr').or.&
                cline%defined('sherr') ) call b%a%introd_alig_err(p%angerr, p%sherr)
            if( p%ctf .eq. 'yes' )       call b%a%introd_ctf_err(p%dferr)
        endif
        if( cline%defined('e1') )then ! rotate input Eulers
            call orientation%new
            call orientation%set_euler([p%e1,p%e2,p%e3])
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
        if( p%zero .eq. 'yes' ) call b%a%zero_shifts
        if( p%discrete .eq. 'yes' )then
            if( cline%defined('ndiscrete') )then
                call b%a%discretize(p%ndiscrete)
            else
                stop 'need ndiscrete to be defined!'
            endif
        endif
        if( cline%defined('nstates') ) call b%a%rnd_states(p%nstates)
        call binwrite_oritab(p%outfile, b%a, [1,b%a%get_noris()])
        call simple_end('**** SIMPLE_ORISOPS NORMAL STOP ****')
    end subroutine exec_orisops

    !> for analyzing SIMPLE orientation/parameter files
    subroutine exec_oristats(self, cline)
        class(oristats_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(build)          :: b
        type(oris)           :: o, osubspace
        type(ori)            :: o_single
        type(params)         :: p
        real                 :: mind, maxd, avgd, sdevd, sumd, vard, scale
        real                 :: mind2, maxd2, avgd2, sdevd2, vard2
        real                 :: popmin, popmax, popmed, popave, popsdev, popvar, frac_populated, szmax
        integer              :: nprojs, iptcl, icls, j 
        integer              :: noris, ncls
        real,    allocatable :: pops(:), tmp(:), clustszs(:), sdevs(:)
        integer, allocatable :: clustering(:)
        logical, allocatable :: ptcl_mask(:)
        integer, parameter   :: hlen=50
        logical              :: err
        p = params(cline)
        call b%build_general_tbox(p, cline, do3d=.false.)
        if( cline%defined('oritab2') )then
            ! Comparison
            if( .not. cline%defined('oritab') ) stop 'need oritab for comparison'
            if( binread_nlines(p%oritab) .ne. binread_nlines(p%oritab2) )then
                stop 'inconsistent number of lines in the two oritabs!'
            endif
            o = oris(p%nptcls)
            call binread_oritab(p%oritab2, o, [1,p%nptcls])
            call b%a%diststat(o, sumd, avgd, sdevd, mind, maxd)
            write(*,'(a,1x,f15.6)') 'SUM OF ANGULAR DISTANCE BETWEEN ORIENTATIONS  :', sumd
            write(*,'(a,1x,f15.6)') 'AVERAGE ANGULAR DISTANCE BETWEEN ORIENTATIONS :', avgd
            write(*,'(a,1x,f15.6)') 'STANDARD DEVIATION OF ANGULAR DISTANCES       :', sdevd
            write(*,'(a,1x,f15.6)') 'MINIMUM ANGULAR DISTANCE                      :', mind
            write(*,'(a,1x,f15.6)') 'MAXIMUM ANGULAR DISTANCE                      :', maxd
        else if( cline%defined('oritab') )then
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
            if( p%classtats .eq. 'yes' )then
                noris = b%a%get_noris()
                ncls  = b%a%get_n('class')
                ! setup weights
                if( p%weights2D.eq.'yes' )then
                    if( noris <= SPECWMINPOP )then
                        call b%a%set_all2single('w', 1.0)
                    else
                        ! frac is one by default in prime2D (no option to set frac)
                        ! so spectral weighting is done over all images
                        call b%a%calc_spectral_weights(1.0)
                    endif
                else
                    ! defaults to unitary weights
                    call b%a%set_all2single('w', 1.0)
                endif
                ! generate class stats
                pops           = b%a%get_pops('class', consider_w=.true.)  !! realloc warning
                popmin         = minval(pops)
                popmax         = maxval(pops)
                popmed         = median_nocopy(pops)
                call moment(pops, popave, popsdev, popvar, err)
                write(*,'(a,1x,f8.2)') 'MINIMUM POPULATION :', popmin
                write(*,'(a,1x,f8.2)') 'MAXIMUM POPULATION :', popmax
                write(*,'(a,1x,f8.2)') 'MEDIAN  POPULATION :', popmed
                write(*,'(a,1x,f8.2)') 'AVERAGE POPULATION :', popave
                write(*,'(a,1x,f8.2)') 'SDEV OF POPULATION :', popsdev
                ! produce a histogram of class populations
                szmax = maxval(pops)
                ! scale to max 50 *:s
                scale = 1.0
                do while( nint(scale*szmax) > hlen )
                    scale = scale - 0.001
                end do
                write(*,'(a)') '>>> HISTOGRAM OF CLASS POPULATIONS'
                do icls=1,ncls
                    write(*,*) nint(pops(icls)),"|",('*', j=1,nint(pops(icls)*scale))  
                end do
            endif
            if( p%projstats .eq. 'yes' )then
                if( .not. cline%defined('nspace') ) stop 'need nspace command line arg to provide projstats'
                noris = b%a%get_noris()
                ! setup weights
                if( p%weights3D.eq.'yes' )then
                    if( noris <= SPECWMINPOP )then
                        call b%a%calc_hard_weights(p%frac)
                    else
                        call b%a%calc_spectral_weights(p%frac)
                    endif
                else
                    call b%a%calc_hard_weights(p%frac)
                endif
                ! generate population stats
                tmp            = b%a%get_pops('proj', consider_w=.true.)!! realloc warning
                nprojs         = size(tmp)
                pops           = pack(tmp, tmp > 0.5)                   !! realloc warning
                frac_populated = real(size(pops))/real(p%nspace)
                popmin         = minval(pops)
                popmax         = maxval(pops)
                popmed         = median_nocopy(pops)
                call moment(pops, popave, popsdev, popvar, err)
                write(*,'(a)') '>>> STATISTICS BEFORE CLUSTERING'
                write(*,'(a,1x,f8.2)') 'FRAC POPULATED DIRECTIONS :', frac_populated
                write(*,'(a,1x,f8.2)') 'MINIMUM POPULATION        :', popmin
                write(*,'(a,1x,f8.2)') 'MAXIMUM POPULATION        :', popmax
                write(*,'(a,1x,f8.2)') 'MEDIAN  POPULATION        :', popmed
                write(*,'(a,1x,f8.2)') 'AVERAGE POPULATION        :', popave
                write(*,'(a,1x,f8.2)') 'SDEV OF POPULATION        :', popsdev
                ! produce a histogram based on clustering into NSPACE_BALANCE even directions
                ! first, generate a mask based on state flag and w
                ptcl_mask = b%a%included(consider_w=.true.)
                allocate(clustering(noris), clustszs(NSPACE_BALANCE))
                call osubspace%new(NSPACE_BALANCE)
                call osubspace%spiral(p%nsym, p%eullims)
                call binwrite_oritab('even_pdirs'//METADATEXT, osubspace, [1,NSPACE_BALANCE])
                do iptcl=1,b%a%get_noris()
                    if( ptcl_mask(iptcl) )then
                        o_single = b%a%get_ori(iptcl)
                        clustering(iptcl) = osubspace%find_closest_proj(o_single)
                    else
                        clustering(iptcl) = 0
                    endif
                end do
                ! determine cluster sizes
                do icls=1,NSPACE_BALANCE
                    clustszs(icls) = real(count(clustering == icls))
                end do
                frac_populated = real(count(clustszs > 0.5))/real(NSPACE_BALANCE)
                popmin         = minval(clustszs)
                popmax         = maxval(clustszs)
                popmed         = median_nocopy(clustszs)
                call moment(clustszs, popave, popsdev, popvar, err)
                write(*,'(a)') '>>> STATISTICS AFTER CLUSTERING'
                write(*,'(a,1x,f8.2)') 'FRAC POPULATED DIRECTIONS :', frac_populated
                write(*,'(a,1x,f8.2)') 'MINIMUM POPULATION        :', popmin
                write(*,'(a,1x,f8.2)') 'MAXIMUM POPULATION        :', popmax
                write(*,'(a,1x,f8.2)') 'MEDIAN  POPULATION        :', popmed
                write(*,'(a,1x,f8.2)') 'AVERAGE POPULATION        :', popave
                write(*,'(a,1x,f8.2)') 'SDEV OF POPULATION        :', popsdev
                ! scale to max 50 *:s
                scale = 1.0
                do while( nint(scale * popmax) > hlen )
                    scale = scale - 0.001
                end do
                write(*,'(a)') '>>> HISTOGRAM OF SUBSPACE POPULATIONS (FROM NORTH TO SOUTH)'
                do icls=1,NSPACE_BALANCE
                    write(*,*) nint(clustszs(icls)),"|",('*', j=1,nint(clustszs(icls)*scale))  
                end do
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
            if( cline%defined('sdev_thres') )then
                sdevs = b%a%get_all('sdev')
                write(*,'(a,1x,i9)') '# particles included:', count(sdevs <= p%sdev_thres)
            endif
        endif
        999 call simple_end('**** SIMPLE_ORISTATS NORMAL STOP ****')
    end subroutine exec_oristats

    !> convert rotation matrix to orientation oris class
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
        call binwrite_oritab(p%outfile, os_out, [1,ndatlines])
        call rotmats%kill
        call simple_end('**** ROTMATS2ORIS NORMAL STOP ****')
    end subroutine exec_rotmats2oris

    !> convert text (.txt) oris doc to binary (.bin)
    subroutine exec_txt2bin( self, cline )
        use simple_oris, only: oris
        class(txt2bin_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)  :: p
        type(oris)    :: os
        integer       :: noris
        p = params(cline)
        noris = nlines(p%oritab)
        call os%new(noris)
        call os%read(p%oritab)
        call binwrite_oritab(p%outfile, os, [1,noris])
    end subroutine exec_txt2bin

    !> convert binary (.bin) oris doc to text (.txt)
    subroutine exec_bin2txt( self, cline )
        use simple_oris, only: oris
        class(bin2txt_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)  :: p
        type(oris)    :: os
        integer       :: noris
        p = params(cline)
        noris = binread_nlines(p%oritab)
        call os%new(noris)
        call binread_oritab(p%oritab, os, [1,noris])
        call os%write(p%outfile)
    end subroutine exec_bin2txt

    subroutine exec_vizoris( self, cline )
        use simple_oris,      only: oris
        use simple_ori,       only: ori
        class(vizoris_commander),  intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(build)           :: b
        type(params)          :: p
        type(ori)             :: o, o_prev
        real,    allocatable  :: euldists(:)
        integer, allocatable  :: pops(:)
        character(len=STDLEN) :: fname, ext
        integer               :: i, n, maxpop, funit, closest,io_stat
        real                  :: radius, maxradius, ang, scale, col, avg_geodist,avg_euldist,geodist
        real                  :: xyz(3), xyz_end(3), xyz_start(3), vec(3)
        p = params(cline)
        call b%build_general_tbox(p, cline, do3d=.false.)
        ! BELOW IS FOR TESTING ONLY
        ! call b%a%spiral
        ! call a%new(2*b%a%get_noris()-1)
        ! do i = 1, b%a%get_noris()
        !     call a%set_ori(i,b%a%get_ori(i))
        ! enddo
        ! do i = b%a%get_noris()+1, 2*b%a%get_noris()-1
        !     call a%set_ori(i,b%a%get_ori(2*b%a%get_noris()-i))
        ! enddo
        ! b%a = a
        n = b%a%get_noris()
        if( .not.cline%defined('fbody') )then
            fname = remove_abspath(trim(adjustl(p%oritab)))
            ext   = trim(fname2ext(fname))
            p%fbody = trim(get_fbody(trim(fname), trim(ext)))
        endif
        if( p%tseries.eq.'no' )then
            ! Discretization of the projection directions
            ! init
            allocate(pops(p%nspace), source=0,stat=alloc_stat)
            allocchk("In commander_oris:: vizoris allocating pops ")
            ang = 3.6 / sqrt(real(p%nsym*p%nspace))
            maxradius = 0.75 * sqrt( (1.-cos(ang))**2. + sin(ang)**2. )
            ! projection direction attribution
            n = b%a%get_noris()
            do i = 1, n
                o = b%a%get_ori(i)
                if(nint(o%get('state')) == 0 )cycle
                call progress(i, n)
                closest = b%e%find_closest_proj(o)
                pops(closest) = pops(closest) + 1          
            enddo        
            maxpop = maxval(pops)
            write(*,'(A,I6)')'>>> NUMBER OF POPULATED PROJECTION DIRECTIONS:', count(pops>0)
            write(*,'(A,I6)')'>>> NUMBER OF EMPTY     PROJECTION DIRECTIONS:', count(pops==0)
            ! output
            fname = trim(p%fbody)//'.bild'
            call fopen(funit, status='REPLACE', action='WRITE', file=trim(fname),iostat=io_stat)
            call fileio_errmsg("simple_commander_oris::exec_vizoris fopen failed "//trim(fname), io_stat)
            ! header
            write(funit,'(A)')".translate 0.0 0.0 0.0"
            write(funit,'(A)')".scale 10"
            write(funit,'(A)')".comment -- unit sphere --"
            write(funit,'(A)')".color 0.8 0.8 0.8"
            write(funit,'(A)')".sphere 0 0 0 1.0"
            write(funit,'(A)')".comment -- planes --"
            write(funit,'(A)')".color 0.3 0.3 0.3"
            write(funit,'(A)')".cylinder -0.02 0 0 0.02 0 0 1.02"
            write(funit,'(A)')".cylinder 0 -0.02 0 0 0.02 0 1.02"
            write(funit,'(A)')".cylinder 0 0 -0.02 0 0 0.02 1.02"
            write(funit,'(A)')".comment -- x-axis --"
            write(funit,'(A)')".color 1 0 0"
            write(funit,'(A)')".cylinder -1.5 0 0 1.5 0 0 0.02"
            write(funit,'(A)')".comment -- y-axis --"
            write(funit,'(A)')".color 0 1 0"
            write(funit,'(A)')".cylinder 0 -1.5 0 0 1.5 0 0.02"
            write(funit,'(A)')".comment -- z-axis --"
            write(funit,'(A)')".color 0 0 1"
            write(funit,'(A)')".cylinder 0 0 -1.5 0 0 1.5 0.02"
            ! body
            write(funit,'(A)')".comment -- projection firections --"
            write(funit,'(A)')".color 0.4 0.4 0.4"
            do i = 1, p%nspace
                if( pops(i) == 0 )cycle
                scale     = real(pops(i)) / real(maxpop)
                xyz_start = b%e%get_normal(i)
                xyz_end   = (1.05 + scale/4.) * xyz_start
                radius    = max(maxradius * scale, 0.002)
                write(funit,'(A,F7.3,F7.3,F7.3,F7.3,F7.3,F7.3,F6.3)')&
                &'.cylinder ', xyz_start, xyz_end, radius
            enddo
            call fclose(funit, errmsg="simple_commander_oris::exec_vizoris closing "//trim(fname))
        else
            ! time series
            ! unit sphere tracking
            fname  = trim(p%fbody)//'_motion.bild'
            radius = 0.02
            call fopen(funit, status='REPLACE', action='WRITE', file=trim(fname), iostat=io_stat)
            call fileio_errmsg("simple_commander_oris::exec_vizoris fopen failed ", io_stat)
            write(funit,'(A)')".translate 0.0 0.0 0.0"
            write(funit,'(A)')".scale 1"
            do i = 1, n
                o   = b%a%get_ori(i)
                xyz = o%get_normal()
                col = real(i-1)/real(n-1)
                write(funit,'(A,F6.2,F6.2)')".color 1.0 ", col, col
                if( i==1 )then
                    write(funit,'(A,F7.3,F7.3,F7.3,A)')".sphere ", xyz, " 0.08"
                else
                    vec = xyz - xyz_start
                    if( sqrt(dot_product(vec, vec)) > 0.01 )then
                        write(funit,'(A,F7.3,F7.3,F7.3,A)')".sphere ", xyz, " 0.02"
                        !write(funit,'(A,F7.3,F7.3,F7.3,F7.3,F7.3,F7.3,F6.3)')&
                        !&'.cylinder ', xyz_start, xyz, radius
                    endif
                endif
                xyz_start = xyz
            enddo
            write(funit,'(A,F7.3,F7.3,F7.3,A)')".sphere ", xyz, " 0.08"
            call fclose(funit, errmsg="simple_commander_oris::exec_vizoris closing "//trim(fname))
            ! distance output
            avg_geodist = 0.
            avg_euldist = 0.
            allocate(euldists(n), stat=alloc_stat)
            fname  = trim(p%fbody)//'_motion.csv'
            call fopen(funit, status='REPLACE', action='WRITE', file=trim(fname), iostat=io_stat)
            call fileio_errmsg("simple_commander_oris::exec_vizoris fopen failed "//trim(fname), io_stat)
            do i = 1, n
                o = b%a%get_ori(i)
                if( i==1 )then
                    ang     = 0.
                    geodist = 0.
                else
                    ang     = rad2deg(o_prev.euldist.o)
                    geodist = o_prev.geod.o
                    call o_prev%mirror2d
                    ang     = min(ang, rad2deg(o_prev.euldist.o))
                    geodist = min(geodist, o_prev.geod.o)
                    avg_euldist = avg_euldist + ang
                    avg_geodist = avg_geodist + geodist
                endif
                euldists(i) = ang
                write(funit,'(I7,A1,F8.3,A1,F8.3)')i, ',', ang, ',', geodist   
                o_prev = o
            enddo
            call fclose(funit, errmsg="simple_commander_oris::exec_vizoris closing "//trim(fname))
            avg_geodist = avg_geodist / real(n-1)
            avg_euldist = avg_euldist / real(n-1)
            write(*,'(A,F8.3)')'>>> AVERAGE EULER    DISTANCE: ',avg_euldist
            write(*,'(A,F8.3)')'>>> AVERAGE GEODESIC DISTANCE: ',avg_geodist
            ! movie output
            ! setting Rprev to south pole as it where chimera opens
            ! call o_prev%new
            ! call o_prev%mirror3d
            ! Rprev = o_prev%get_mat()
            ! fname = trim(p%fbody)//'_movie.cmd'
            ! call fopen(funit, status='REPLACE', action='WRITE', file=trim(fname), iostat=io_stat)
            ! call fileio_errmsg("simple_commander_oris::exec_vizoris fopen failed "//trim(fname), io_stat)
            ! do i = 1, n
            !     o  = b%a%get_ori(i)
            !     Ri = o%get_mat()
            !     R  = matmul(Ri,transpose(Rprev))
            !     call rotmat2axis(R, axis)
            !     if( abs(euldists(i)) > 0.01 )then
            !         if(i==1)then
            !             euldists(1) = rad2deg(o.euldist.o_prev)
            !             write(funit,'(A,A,A1,A,A1,A,A1,F8.2,A2)')&
            !             &'roll ', trim(real2str(axis(1))),',', trim(real2str(axis(2))),',', trim(real2str(axis(3))),' ',&
            !             &euldists(i), ' 1'
            !         else
            !             write(funit,'(A,A,A1,A,A1,A,A1,F8.2,A)')&
            !             &'roll ', trim(real2str(axis(1))),',', trim(real2str(axis(2))),',', trim(real2str(axis(3))),' ',&
            !             &euldists(i), ' 2; wait 2'
            !         endif
            !         o_prev = o
            !         Rprev  = Ri
            !     endif
            ! enddo
            ! call fclose(funit, errmsg="simple_commander_oris::exec_vizoris closing "//trim(fname))
        endif
        call simple_end('**** VIZORIS NORMAL STOP ****')
    end subroutine exec_vizoris

end module simple_commander_oris
