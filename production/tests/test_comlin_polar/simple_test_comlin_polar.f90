program simple_test_comlin_polar
include 'simple_lib.f08'
use simple_cmdline,          only: cmdline
use simple_parameters,       only: parameters
use simple_image,            only: image
use simple_projector,        only: projector
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_polarizer,        only: polarizer
use simple_sym,              only: sym
implicit none

type polar_fmap
    logical :: legit
    integer :: targ_irot
    real    :: targ_w
    integer :: self_irot
    real    :: self_w
end type polar_fmap
character(len=*), parameter   :: POLAR_COMLIN = 'polar_comlin.bin'
integer,          parameter   :: NPLANES = 300, ORI_IND1 = 10, ORI_IND2 = 55
character(len=:), allocatable :: cmd
type(image),      allocatable :: pad_fplanes(:)
complex,          allocatable :: cmat(:,:), pfts(:,:,:), even_pfts(:,:,:), line1(:), line2(:)
complex,          pointer     :: ref_ptrs_even(:,:,:), ref_ptrs_odd(:,:,:)
type(parameters)              :: p
type(polarft_corrcalc)        :: pftcc
type(polarizer)               :: img_polarizer
type(cmdline)                 :: cline
type(image)                   :: noise, ptcl, ptcl_pad, fplanes(NPLANES), vol, img, fplane_polar
type(oris)                    :: spiral
type(ori)                     :: o1, o2
type(projector)               :: vol_pad
type(polar_fmap)              :: pcomlines(NPLANES,NPLANES)
type(sym)                     :: se
integer :: ifoo, rc, i, lims(3,2), ithr, box, errflg, kfromto(2), j,k
real    :: ave, sdev, maxv, minv, ang, ang1, d, dang
logical :: mrc_exists, lineinterp
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)') 'ERROR! Usage: simple_test_comlin_polar smpd=xx nthr=yy vol1=volume.mrc mskdiam=zz'
    write(logfhandle,'(a)') 'Example: https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
    write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
    inquire(file="1JYX.mrc", exist=mrc_exists)
    if( .not. mrc_exists )then
        write(*, *) 'Downloading the example dataset...'
        cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
        call execute_command_line(cmd, exitstat=rc)
        write(*, *) 'Converting .pdb to .mrc...'
        cmd = 'e2pdb2mrc.py 1JYX.pdb 1JYX.mrc'
        call execute_command_line(cmd, exitstat=rc)
        cmd = 'rm 1JYX.pdb'
        call execute_command_line(cmd, exitstat=rc)
    endif
    call cline%set('smpd'   , 1.)
    call cline%set('nthr'   , 16.)
    call cline%set('vol1'   , '1JYX.mrc')
    call cline%set('mskdiam', 180.)
    call cline%set('lp'   ,   3.)
else
    call cline%parse_oldschool
endif
call cline%checkvar('smpd',    1)
call cline%checkvar('nthr',    2)
call cline%checkvar('vol1',    3)
call cline%checkvar('mskdiam', 4)
call cline%checkvar('lp',      5)
call cline%check
call p%new(cline)
call find_ldim_nptcls(p%vols(1), p%ldim, ifoo)
call vol%new(p%ldim, p%smpd)
call vol%read(p%vols(1))
call vol%stats('foreground', ave, sdev, maxv, minv)
call spiral%new(NPLANES, is_ptcl=.false.)
se = sym('c1')
call se%build_refspiral(spiral)
call ptcl%new(    [p%box,   p%box,   1],       p%smpd)
call noise%new(   [p%box,   p%box ,  p%box],   p%smpd)
call img%new(     [p%box,   p%box,   1],       p%smpd)
call vol_pad%new( [p%boxpd, p%boxpd, p%boxpd], p%smpd)
call ptcl_pad%new([p%boxpd, p%boxpd, 1],       p%smpd)
! add noise in a small center region of the vol
call noise%gauran(0., 5. * sdev)
call noise%mask(1.5 * p%msk, 'soft')
call vol%add(noise)
call vol%pad(vol_pad)
call vol_pad%fft
call vol_pad%expand_cmat(p%alpha)
call spiral%get_ori(ORI_IND1, o1)
call vol_pad%fproject(o1,ptcl_pad)
call ptcl_pad%ifft
call ptcl_pad%clip(ptcl)
call ptcl%write('reproj_com_reprojcom.mrc', 1)
lims = vol%loop_lims(3)
allocate(pad_fplanes(p%nthr))
do ithr = 1, p%nthr
    call pad_fplanes(ithr)%new([p%boxpd, p%boxpd, 1], p%smpd)
enddo
!$omp parallel do default(shared) private(i,ithr,o2)&
!$omp proc_bind(close) schedule(static)
do i = 1, spiral%get_noris()
    ithr = omp_get_thread_num() + 1
    ! fplanes are used for optimization below
    call spiral%get_ori(i, o2)
    call fplanes(i)%new([p%box, p%box, 1], p%smpd)
    call vol_pad%fproject(o2,pad_fplanes(ithr))
    call pad_fplanes(ithr)%ifft
    call pad_fplanes(ithr)%clip(fplanes(i))
enddo
!$omp end parallel do
do i = 1, spiral%get_noris()
    call fplanes(i)%write('fplanes.mrc', i)
    call fplanes(i)%fft
enddo
! polar stuffs
kfromto = p%kfromto
call img_polarizer%new([p%box,p%box,1],p%smpd, wthreads=.false.)
call pftcc%new(NPLANES, [1,NPLANES], kfromto)
call img_polarizer%init_polarizer(pftcc, p%alpha)
call fplane_polar%new([p%box,p%box,1],1.0)
do i = 1, spiral%get_noris()
    call img_polarizer%polarize(pftcc, fplanes(i), i, isptcl=.false., iseven=.true.)
    call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
    call fplane_polar%zero_and_flag_ft
    call fplane_polar%set_cmat(cmat)
    call fplane_polar%shift_phorig()
    call fplane_polar%ifft
    call fplane_polar%write('fplanes_polar.mrc', i)
enddo
! Testing polar line generation
call pftcc%get_refs_ptr(ref_ptrs_even, ref_ptrs_odd)
! Test that reproduces the strategy in gen_polar_comlins
lineinterp = .true.
allocate(line1(kfromto(1):kfromto(2)), line2(kfromto(1):kfromto(2)))
dang = pftcc%get_rot(2)-pftcc%get_rot(1)
do i = 1,pftcc%get_nrots()
    ang1 = pftcc%get_rot(i)     ! angle in pftcc
    ang  = ang1 + ran3()*dang   ! random angle between two lines
    ! line1: line at ang from intepolated pftcc lines
    d    = ang-ang1
    j    = i - pftcc%get_pftsz()
    if( i < pftcc%get_pftsz() )then
        line1 = (1.0-d)*ref_ptrs_even(i,:,1) + d*ref_ptrs_even(i+1,:,1)
    elseif( i == pftcc%get_pftsz() )then
        line1 = (1.0-d)*ref_ptrs_even(i,:,1) + d*conjg(ref_ptrs_even(j+1,:,1))
    elseif( i < pftcc%get_nrots() )then
        line1 = (1.0-d)*conjg(ref_ptrs_even(j,:,1)) + d*conjg(ref_ptrs_even(j+1,:,1))
    else
        line1 = (1.0-d)*conjg(ref_ptrs_even(j,:,1))  + d*ref_ptrs_even(1,:,1)
    endif
    ! line2: line at ang from image pixel interpolation
    do k = kfromto(1), kfromto(2)
        line2(k) = img_polarizer%extract_pixel(ang,k,fplanes(1))
    enddo
    if( sum(csq_fast(line1-line2)) > 1.e-7 )then
        print *,'Line intepolation test failed for line: '
        print *,i,ang1,ang,sum(csq_fast(line1-line2))
        print *,line1
        print *,line2
        lineinterp = .false.
    endif
enddo
if( lineinterp )then
    print *,' TEST FOR INTERPOLATION STRATEGY PASSED'
else
    print *,' TEST FOR INTERPOLATION STRATEGY FAILED'
endif
allocate(pfts(pftcc%get_pftsz(), kfromto(1):kfromto(2), NPLANES), source=cmplx(0.,0.))
allocate(even_pfts(pftcc%get_pftsz(), kfromto(1):kfromto(2), NPLANES), source=ref_ptrs_even) ! backup
call gen_polar_comlins(pftcc, spiral, pcomlines)
call polar_comlin_pfts(pcomlines, ref_ptrs_even, pfts)
i = ORI_IND1
call pftcc%set_ref_pft(i, ref_ptrs_even(:,:,i), iseven=.true.) ! not necessay as ref_ptrs_even => pftcc%refs_pfts_even
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%zero_and_flag_ft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 1)
call pftcc%set_ref_pft(i, pfts(:,:,i), iseven=.true.) ! this sets refs_pfts_even to pfts
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%zero_and_flag_ft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 2)
where( sqrt(real(pfts(:,:,i)*conjg(pfts(:,:,i)))) > TINY ) even_pfts(:,:,i) = (even_pfts(:,:,i) + pfts(:,:,i))/2.
call pftcc%set_ref_pft(i, even_pfts(:,:,i), iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%zero_and_flag_ft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 3)
i = ORI_IND2
call pftcc%set_ref_pft(i, ref_ptrs_even(:,:,i), iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%zero_and_flag_ft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 4)
call pftcc%set_ref_pft(i, pfts(:,:,i), iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%zero_and_flag_ft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 5)
where( sqrt(real(pfts(:,:,i)*conjg(pfts(:,:,i)))) > TINY ) even_pfts(:,:,i) = (even_pfts(:,:,i) + pfts(:,:,i))/2.
call pftcc%set_ref_pft(i, even_pfts(:,:,i), iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%zero_and_flag_ft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 6)

contains

    subroutine get_polar_coord( pftcc, xy, rot, k )
        class(polarft_corrcalc), intent(in)  :: pftcc
        real,                    intent(in)  :: xy(2)
        real,                    intent(out) :: rot, k
        real :: tan_inv
        k       = sqrt(sum(xy**2))
        tan_inv = atan(xy(2), xy(1)) * 2 + PI
        rot     = 1. + tan_inv * real(pftcc%get_pftsz()) / TWOPI
    end subroutine get_polar_coord

    subroutine gen_polar_comlins( pftcc, ref_space, pcomlines)
        class(polarft_corrcalc), intent(in)    :: pftcc
        type(oris),              intent(in)    :: ref_space
        type(polar_fmap),        intent(inout) :: pcomlines(:,:)
        logical, parameter :: l_kb = .true. ! K-B or linear interpolation
        integer :: iref, jref, irot, kind, irot_l, nrefs, nrots
        real    :: loc1_3D(3), loc2_3D(3), denom, a1, a2, b1, b2, line2D(3), irot_real, k_real, w, hk1(2), hk2(2),&
                  &rotmat(3,3),invmats(3,3,pftcc%get_nrefs()), loc1s(3,pftcc%get_nrefs()), loc2s(3,pftcc%get_nrefs()), line3D(3),&
                  &euls(3,pftcc%get_nrefs()), euldist
        logical :: l_exclude_cone
        nrefs = pftcc%get_nrefs()
        nrots = pftcc%get_nrots()
        ! randomly chosing two sets of (irot, kind) to generate the polar common lines
        irot  = 5
        kind  = p%kfromto(1) + 5
        if( kind >= p%kfromto(2) ) kind = p%kfromto(1) ! using kfromto(1) as kind in the first set of (irot,kind)
        hk1   = pftcc%get_coord(irot,kind)
        irot  = 16
        kind  = p%kfromto(1) + 15
        if( kind >= p%kfromto(2) ) kind = p%kfromto(2) ! using kfromto(2) as kind in the second set of (irot,kind)
        hk2   = pftcc%get_coord(irot,kind)
        ! angular threshold
        l_exclude_cone = trim(p%linethres).eq.'yes'
        ! caching rotation matrices and their corresponding inverse matrices
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,rotmat)
        do iref = 1, nrefs
            rotmat            = ref_space%get_mat(iref)
            invmats(:,:,iref) = transpose(rotmat)
            loc1s(:,iref)     = matmul([hk1(1), hk1(2), 0.], rotmat)
            loc2s(:,iref)     = matmul([hk2(1), hk2(2), 0.], rotmat)
            euls( :,iref)     = ref_space%get_euler(iref)
        enddo
        !$omp end parallel do
        euls(3,:) = 0.
        ! constructing polar common lines
        pcomlines%legit = .false.
        !$omp parallel do default(shared) proc_bind(close) schedule(static)&
        !$omp private(iref,loc1_3D,loc2_3D,denom,a1,b1,jref,euldist,a2,b2,line3D,line2D,irot_real,k_real,irot_l,w)
        do iref = 1, nrefs
            loc1_3D = loc1s(:,iref)
            loc2_3D = loc2s(:,iref)
            denom   = (loc1_3D(1) * loc2_3D(2) - loc1_3D(2) * loc2_3D(1))
            a1      = (loc1_3D(3) * loc2_3D(2) - loc1_3D(2) * loc2_3D(3)) / denom
            b1      = (loc1_3D(1) * loc2_3D(3) - loc1_3D(3) * loc2_3D(1)) / denom
            do jref = 1, nrefs
                if( jref == iref )cycle
                if( l_exclude_cone )then
                    euldist = rad2deg(euler_dist(euls(:,iref),euls(:,jref)))
                    if( euldist < p%athres )cycle
                endif
                ! getting the 3D common line
                loc1_3D     = loc1s(:,jref)
                loc2_3D     = loc2s(:,jref)
                denom       = (loc1_3D(1) * loc2_3D(2) - loc1_3D(2) * loc2_3D(1))
                if( abs(denom) < TINY )cycle
                a2          = (loc1_3D(3) * loc2_3D(2) - loc1_3D(2) * loc2_3D(3)) / denom
                b2          = (loc1_3D(1) * loc2_3D(3) - loc1_3D(3) * loc2_3D(1)) / denom
                if( abs(b1-b2) < TINY )cycle
                line3D(1:2) = [1., -(a1-a2)/(b1-b2)]
                line3D(3)   = a2*line3D(1) + b2*line3D(2)
                ! projecting the 3D common line to a polar line on the jref-th reference
                line2D      = matmul(line3D, invmats(:,:,jref))
                call get_polar_coord(pftcc, line2D(1:2), irot_real, k_real)
                ! caching the indices irot_j and irot_j+1 and the corresponding linear weight
                if( l_kb )then
                    ! nearest grid point
                    irot_l = nint(irot_real)
                else
                    ! leftmost grid point
                    irot_l = floor(irot_real)
                endif
                w = irot_real - real(irot_l)
                pcomlines(jref,iref)%targ_irot = irot_l
                pcomlines(jref,iref)%targ_w    = w
                ! projecting the 3D common line to a polar line on the iref-th reference
                line2D = matmul(line3D, invmats(:,:,iref))
                call get_polar_coord(pftcc, line2D(1:2), irot_real, k_real)
                ! caching the indices irot_i and irot_i+1 and the corresponding linear weight
                if( l_kb )then
                    ! nearest grid point
                    irot_l = nint(irot_real)
                else
                    ! leftmost grid point
                    irot_l = floor(irot_real)
                endif
                w = irot_real - real(irot_l)
                pcomlines(jref,iref)%self_irot = irot_l
                pcomlines(jref,iref)%self_w    = w
                pcomlines(jref,iref)%legit     = .true.
            enddo
        enddo
        !$omp end parallel do
    end subroutine gen_polar_comlins

    subroutine gen_polar_comlinsm( pftcc, ref_space, pcomlines)
        class(polarft_corrcalc), intent(in)    :: pftcc
        type(oris),              intent(in)    :: ref_space
        type(polar_fmap),        intent(inout) :: pcomlines(:,:)
        logical, parameter :: l_kb = .true. ! K-B or linear interpolation
        real    :: eulers(3), euldist, R(3,3,pftcc%get_nrefs()), Rj(3,3), tRi(3,3)
        real    :: psii, psij, drot, d
        integer :: iref, jref, irot, nrefs, nrots, imirr
        logical :: l_exclude_cone
        if( .not.ref_space%isthere('mirr') )then
            print *, 'Mirror index missing in reference search space'
            stop
        endif
        nrefs = pftcc%get_nrefs()
        nrots = pftcc%get_nrots()
        drot  = 360.0 / real(nrots)
        ! angular threshold
        l_exclude_cone = trim(p%linethres).eq.'yes'
        !$omp parallel default(shared) proc_bind(close)&
        !$omp private(iref,jref,imirr,tRi,euldist,Rj,eulers,irot,d,psii,psij)
        ! caching rotation matrices
        !$omp do schedule(static)
        do iref = 1, nrefs
            R(:,:,iref) = ref_space%get_mat(iref)
        enddo
        !$omp end do
        ! constructing polar common lines
        !$omp do schedule(static)
        do iref = 1,nrefs
            tRi   = transpose(R(:,:,iref))
            imirr = ref_space%get_int(iref,'mirr')
            do jref = 1,nrefs
                pcomlines(jref,iref)%legit  = .false.
                if( jref == iref  ) cycle   ! self   exclusion
                if( jref == imirr ) cycle   ! mirror exclusion
                ! Cone exclusion
                if( l_exclude_cone )then
                    euldist = rad2deg(ref_space%euldist(iref,jref))
                    if( euldist < p%athres )cycle
                endif
                ! Rotation of both planes by transpose of Ri such plane i has euler triplet (0,0,0)
                Rj     = matmul(R(:,:,jref), tRi)
                eulers = m2euler(Rj)
                ! in plane rotation index of jref plane intersecting iref
                psij = 360.0 - eulers(3)
                irot = pftcc%get_roind(psij)
                d    = psij - pftcc%get_rot(irot)
                if( d > drot ) d = d - 360.0
                if( l_kb )then
                    ! nearest point with get_roind routine
                else
                    ! leftmost point
                    if( d < 0. )then
                        irot = irot - 1
                        d    = d + drot
                        if( irot < 1 ) irot = irot + nrots
                    endif
                endif
                pcomlines(jref,iref)%targ_irot = irot
                pcomlines(jref,iref)%targ_w    = d / drot
                ! in plane rotation index of iref plane (ccw)
                psii = eulers(1)
                irot = pftcc%get_roind(psii)
                d    = psii - pftcc%get_rot(irot)
                if( d > drot ) d = d - 360.0
                if( l_kb )then
                    ! nearest point with get_roind routine
                else
                    ! leftmost point
                    if( d < 0. )then
                        irot = irot - 1
                        d    = d + drot
                        if( irot < 1 ) irot = irot + nrots
                    endif
                endif
                pcomlines(jref,iref)%self_irot = irot
                pcomlines(jref,iref)%self_w    = d / drot
                ! green lighting pair
                pcomlines(jref,iref)%legit = .true.
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine gen_polar_comlinsm

    subroutine polar_comlin_pfts( pcomlines, pfts_in, pfts)
        type(polar_fmap), intent(in)    :: pcomlines(:,:)
        complex,          intent(in)    :: pfts_in(:,:,:)
        complex,          intent(inout) :: pfts(:,:,:)
        complex :: pft_line(p%kfromto(1):p%kfromto(2))
        integer :: iref, jref, irot_l, irot_r, nrefs, pftsz
        real    :: w
        pftsz = size(pfts_in, 1)
        nrefs = size(pcomlines,1)
        pfts  = CMPLX_ZERO
        !$omp parallel do default(shared) private(iref,jref,irot_l,irot_r,w,pft_line)&
        !$omp proc_bind(close) schedule(static)
        do iref = 1, nrefs
            do jref = 1, nrefs
                if( .not. pcomlines(jref,iref)%legit )cycle
                ! compute the interpolated polar common line, between irot_j and irot_j+1
                irot_l = pcomlines(jref,iref)%targ_irot
                irot_r = irot_l + 1
                w      = pcomlines(jref,iref)%targ_w
                if( irot_l < 1 )then
                    pft_line = (1.-w) * conjg(pfts_in(irot_l+pftsz,:,jref))
                elseif( irot_l > pftsz )then
                    pft_line = (1.-w) * conjg(pfts_in(irot_l-pftsz,:,jref))
                else
                    pft_line = (1.-w) * pfts_in(irot_l,:,jref)
                endif
                if( irot_r < 1 )then
                    pft_line = pft_line + w * conjg(pfts_in(irot_r+pftsz,:,jref))
                elseif( irot_r > pftsz )then
                    pft_line = pft_line + w * conjg(pfts_in(irot_r-pftsz,:,jref))
                else
                    pft_line = pft_line + w * pfts_in(irot_r,:,jref)
                endif
                ! extrapolate the interpolated polar common line to irot_i and irot_i+1 of iref-th reference
                irot_l   = pcomlines(jref,iref)%self_irot
                irot_r   = irot_l + 1
                w        = pcomlines(jref,iref)%self_w
                if( irot_l < 1 )then
                    irot_l              = irot_l + pftsz
                    pfts(irot_l,:,iref) = pfts(irot_l,:,iref) + (1.-w) * conjg(pft_line)
                elseif( irot_l > pftsz )then
                    irot_l              = irot_l - pftsz
                    pfts(irot_l,:,iref) = pfts(irot_l,:,iref) + (1.-w) * conjg(pft_line)
                else
                    pfts(irot_l,:,iref) = pfts(irot_l,:,iref) + (1.-w) * pft_line
                endif
                if( irot_r < 1 )then
                    irot_r              = irot_r + pftsz
                    pfts(irot_r,:,iref) = pfts(irot_r,:,iref) + w * conjg(pft_line)
                elseif( irot_r > pftsz )then
                    irot_r              = irot_r - pftsz
                    pfts(irot_r,:,iref) = pfts(irot_r,:,iref) + w * conjg(pft_line)
                else
                    pfts(irot_r,:,iref) = pfts(irot_r,:,iref) + w * pft_line
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine polar_comlin_pfts

    subroutine read_write_comlin( pcomlines, pftcc, eulspace )
        type(polar_fmap), allocatable, intent(inout) :: pcomlines(:,:)
        class(polarft_corrcalc),       intent(in)    :: pftcc
        type(oris),                    intent(in)    :: eulspace
        integer :: io_stat, funit
        if( .not. allocated(pcomlines) ) allocate(pcomlines(p%nspace,p%nspace))
        if( file_exists(trim(POLAR_COMLIN)) )then
            call fopen(funit,trim(POLAR_COMLIN),access='STREAM',action='READ',status='OLD', iostat=io_stat)
            read(unit=funit,pos=1) pcomlines
        else
            call gen_polar_comlinsm(pftcc, eulspace, pcomlines)
            call fopen(funit,trim(POLAR_COMLIN),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
            write(unit=funit,pos=1) pcomlines
        endif
        call fclose(funit)
    end subroutine read_write_comlin

end program simple_test_comlin_polar