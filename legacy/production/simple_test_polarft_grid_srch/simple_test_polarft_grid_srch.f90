program simple_test_polarft_grid_srch
use simple_build,      only: build
use simple_params,     only: params
use simple_polarft,    only: polarft
use simple_align_pair, only: align_pair
use simple_cmdline     ! singleton
use simple_jiffys      ! singleton
use simple_defs        ! singleton
use simple_rnd, only:  irnd_uni, ran3
use simple_polarft_shsrch
implicit none
type(params)          :: p
type(build)           :: b
type(polarft), target :: pft, pft_ref
integer               :: i, j, k, x, y, r, trs, nerrs, loc(1), cnt, cnt2
real                  :: shvec(2), cxy(3), lims(2,2), rot, ang, dist, adist
integer, parameter    :: NROUNDS=5
real, allocatable     :: ransh(:,:), corrs(:), shifts(:,:), rots(:)
if( command_argument_count() < 4 )then
    write(*,'(a)', advance='no') 'SIMPLE_TEST_POLARFT_GRID_SRCH stk=stack.spi smpd=<sampling distance(in A)>'
    write(*,'(a)') ' msk=<mask radius(in pixels)> trs=<origin shift(in pixels)'
    stop
endif
call parse_cmdline
call cmdcheckvar('stk',  1)
call cmdcheckvar('smpd', 2)
call cmdcheckvar('msk',  3)
call cmdcheckvar('trs',  4)
call cmdcheck
p = params()                 ! parameters generated
call b%build_general_tbox(p) ! general objects built
p%fromp = 1+5
p%top   = p%nptcls+5
allocate( ransh(p%nptcls,2), corrs(p%nptcls), shifts(p%nptcls,2), rots(p%nptcls) )
call b%pftcc%new(p%nptcls, [p%fromp,p%top], [p%box,p%box,1], [2,20], p%ring2)
p%kfromto = [2,20]
trs = nint(p%trs)
call pft%new(p%kfromto, p%ring2, ptcl=.true.)
call pft_ref%new(p%kfromto, p%ring2, ptcl=.false.)
lims(:,2) = p%trs
lims(:,1) = -p%trs
call polarft_shsrch_init(pft_ref, pft, lims)
call polarft_shsrch_set_rot(1)
do j=1,NROUNDS
    shvec = [ran3()*p%trs*2.-p%trs,ran3()*p%trs*2.-p%trs]
    do i=1,p%nptcls
        call b%img%read(p%stk, i)
        call b%img%fwd_ft
        call b%proj%img2polarft(b%img, p%msk, pft_ref)
        call b%img%shift(shvec(1), shvec(2),imgout=b%img_copy)
        call b%proj%img2polarft(b%img_copy, p%msk, pft)
        cxy = polarft_shsrch_minimize()
        if( sqrt(sum((shvec-cxy(2:3))**2.)) > 0.12 )then
            stop 'shift-only search failed'
        endif
    end do
end do
deallocate(corrs)
allocate(corrs(pft%get_nrots()))
adist = 0.
do j=1,NROUNDS
    shvec = [ran3()*p%trs*2.-p%trs,ran3()*p%trs*2.-p%trs]
    do i=1,p%nptcls
        call b%img%read(p%stk, i)
        call b%img%fwd_ft
        call b%proj%img2polarft(b%img, p%msk, pft_ref)
        ang = ran3()*360.
        call b%img%rtsq(ang, 0., 0.)
        call b%img%shift(shvec(1), shvec(2),imgout=b%img_copy)
        call b%proj%img2polarft(b%img_copy, p%msk, pft)
        call pft_ref%gencorrs(pft, corrs)
        loc = maxloc(corrs)
        rot = pft%get_rot(loc(1))
        call polarft_shsrch_set_rot(loc(1))
        cxy = polarft_shsrch_minimize()
        dist = sqrt(sum((shvec-cxy(2:3))**2.))
        adist = adist+dist
    end do
end do
adist = adist/real(NROUNDS*p%nptcls)
print *, 'average distance between correct and identified shift vector (in pixels): ', adist
deallocate(corrs)
allocate(corrs(p%nptcls))
nerrs = 0
cnt = 0
do j=1,NROUNDS
    ! make random shift array
    do i=1,p%nptcls
        x = (irnd_uni(trs+1)-1)*2-trs
        y = (irnd_uni(trs+1)-1)*2-trs
        ransh(i,1) = real(x)
        ransh(i,2) = real(y)
    end do
    ! prepare data structure
    cnt2 = 0
    do i=p%fromp,p%top
        cnt2 = cnt2+1
        call b%img%read(p%stk, cnt2)
        call b%img%fwd_ft
        call b%proj%img2polarft(cnt2, b%img, p%msk, b%pftcc, isptcl=.false.) ! this is the reference
        call b%img%shift(ransh(cnt2,1), ransh(cnt2,2), imgout=b%img_copy)
        call b%proj%img2polarft(i, b%img_copy, p%msk, b%pftcc, isptcl=.true.) ! this is the particle (shifted)
    end do
    ! do the grid search
    cnt2 = 0
    do i=p%fromp,p%top
        cnt2 = cnt2+1
        cnt  = cnt+1
        call progress(cnt, NROUNDS*p%nptcls)
!         print *, 'search for particle: ', i
        do k=1,p%nptcls
            call b%pftcc%grid_srch(k, i, real(trs), r, shifts(k,:), corrs(k))
            rots(k) = b%pftcc%get_rot(r)
!             print *, 'ref: ', k, 'corr: ', corrs(k), 'rot: ', rots(k)
        end do
        loc = maxloc(corrs)
!         print *, 'best matching reference: ', loc(1), 'corr: ', corrs(loc(1))
        if( loc(1) /= i-5 )then
            print *, 'identified the wrong reference: ', loc(1)
            stop
        endif
        if( rots(loc(1)) > 1e-6 )then
            print *, 'in-plane rotation is not zero: ', rots(loc(1))
            stop
        endif
!         print *, 'correct origin shift: ', ransh(loc(1),:), 'identified origin shift: ', shifts(loc(1),:)
        if( sqrt(sum((ransh(loc(1),:)-shifts(loc(1),:))**2.)) > 1e-6 )then
            nerrs = nerrs+1
        endif
    end do
end do
print *, 'total number of errors (%)', 100.*(real(nerrs)/real((trs+1)**2*p%nptcls))
call simple_end('**** SIMPLE_TEST_POLARFT_GRID_SRCH NORMAL STOP ****')
end program simple_test_polarft_grid_srch
    