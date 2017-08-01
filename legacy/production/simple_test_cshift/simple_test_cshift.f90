program simple_test_cshift
use simple_image,  only: image
use simple_math,   only: round2even
use simple_rnd,    only: ran3
use simple_jiffys, only: progress
use simple_timing
implicit none

! IDEA: we can safely double the particle's nrots but we better dynamically reshape the array in the first dimension

! type pft_single
!     complex, allocatable :: cmat(:,:)
! end type

! type(pft_single), allocatable :: pfts_ptcls_singles(:)
integer, allocatable :: fake_pft(:,:,:), ptcl_order(:)
integer, allocatable :: roshmat(:,:), ptshmat(:,:)
complex, allocatable :: pfts_ptcls(:,:,:), pfts_ptcls_transf(:,:,:)
integer, parameter   :: KFAKE=4, NROTS_FAKE=3, NREFS_FAKE=9, fromp=1
integer, parameter   :: box=240, nspace=100, ring2=110
real, parameter      :: smpd=1.77, lp=8
type(image)          :: img
integer :: irot, jrot, iref, jref, cnt, cnt2, lim(2), iptcl, jptcl 
integer :: top, nrots, kfromto(2), refsz, ptclsz, i, ptcllim(2), rotlim(2)
integer :: roshiftmat(NREFS_FAKE,KFAKE), ptshiftmat(NROTS_FAKE,KFAKE)

call start_Alltimers_cpu()
top = NREFS_FAKE
allocate( fake_pft(NREFS_FAKE,NROTS_FAKE,KFAKE) )
fake_pft = 0
cnt = 0
! fill-up the fake pft
do iref=1,NREFS_FAKE
    cnt = cnt+10
    cnt2 = 0
    do irot=1,NROTS_FAKE
        cnt2 = cnt2+1
        fake_pft(iref,irot,:) = cnt+cnt2
    end do
end do
! call print_fake_pft
! write (*,*) '>>> TESTING THE MAT PT CHANGER'
! ptshiftmat=1
! fake_pft = cshift(fake_pft, shift=ptshiftmat, dim=1)
! fake_pft = cshift(fake_pft, shift=-ptshiftmat, dim=1)
! call print_fake_pft
! write (*,*) '>>> TESTING THE MAT RO CHANGER'
! roshiftmat=1
! fake_pft = cshift(fake_pft, shift=roshiftmat, dim=2)
! fake_pft = cshift(fake_pft, shift=-roshiftmat, dim=2)
! call print_fake_pft

! stop

! write (*,*) '>>> TESTING THE PTCL CHANGER'
! call print_fake_pft
! do iptcl=2,NREFS_FAKE
!     fake_pft = cshift(fake_pft, shift=1, dim=1)
!     call print_fake_pft
! end do
! ! one last call to get the original order back
! fake_pft = cshift(fake_pft, 1, 1)
! write (*,*) '>>> ORIGINAL PTCL ORDER'
! call print_fake_pft
!
! write (*,*) '>>> TESTING THE ROTATION CHANGER'
! call print_fake_pft
! rosh = 1
! fake_pft = cshift(fake_pft, shift=rosh, dim=2)



! do irot=2,NROTS_FAKE
!     fake_pft = cshift(fake_pft, 1, 2)
!     call print_fake_pft
! end do
! ! one last call to get the original order back
! fake_pft = cshift(fake_pft, 1, 2)
! write (*,*) '>>> ORIGINAL ROTATION ORDER'
! call print_fake_pft

! simulate a real cshift round
nrots  = round2even(twopi*real(ring2))
top    = NSPACE
refsz  = nrots/2
ptclsz = 2*nrots
! ptclsz = nrots
call img%new([box,box,1], smpd)
kfromto(1)  = 2
kfromto(2)  = img%get_find(lp)
allocate(pfts_ptcls(2*nspace,ptclsz,kfromto(1):kfromto(2)), pfts_ptcls_transf(nspace,refsz,kfromto(1):kfromto(2)))
! roshmat = 1
! ptshmat = 1
! call start_timer_cpu("cshift_ver1")
! call gettimeofday_c(st_r)
! do irot=1,nrots
!     call progress(irot, nrots)
!     do iref=1,nspace
!         if( iref /= 1 ) ptcl_order = cshift(ptcl_order, 1, 1)
!         do i=1,nspace
!             iptcl = ptcl_order(i)
!             if( irot /= 1 ) pfts_ptcls_singles(iptcl)%cmat = cshift(pfts_ptcls_singles(iptcl)%cmat, 1, 1)
!             pfts_ptcls_transf(i,:,:) = pfts_ptcls_singles(iptcl)%cmat(1:refsz,:)
!         end do
!         ! calculate hadamard product here
!     end do

! end do
! call gettimeofday_c(et_r)
! call elapsed_time_c(st_r,et_r,elps_r)
! write(*,'(x,a,x,f20.8)') "the elapse time for the tesla_r(s)= ",elps_r
! call stop_timer_cpu("cshift_ver1")



! call start_timer_cpu("cshift_ver2")
! do iptcl=1,nspace
!     call progress(iptcl, nspace)
!     if( iptcl /= 1 ) pfts_ptcls = cshift(pfts_ptcls, shift=ptshmat, dim=1)
!     do irot=1,nrots
!         if( irot /= 1 ) pfts_ptcls = cshift(pfts_ptcls, shift=roshmat, dim=2)
!         ! calculate hadamard product here
!     end do
!     ! one last call to get the original order back
!     pfts_ptcls = cshift(pfts_ptcls, shift=roshmat, dim=2)
! end do
! ! one last call to get the original order back
! pfts_ptcls = cshift(pfts_ptcls, shift=ptshmat, dim=1)
! call stop_timer_cpu("cshift_ver2")
! call start_timer_cpu("cshift_ver3")
! do irot=1,nrots
!     call progress(irot, nrots)
!     if( irot /= 1 ) pfts_ptcls = cshift(pfts_ptcls, 1, 2)
!     do iptcl=1,nspace
!         if( iptcl /= 1 ) pfts_ptcls = cshift(pfts_ptcls, 1, 1)
!         ! calculate hadamard product here
!     end do
!     ! one last call to get the original order back
!     pfts_ptcls = cshift(pfts_ptcls, 1, 1)
! end do
! ! one last call to get the original order back
! pfts_ptcls = cshift(pfts_ptcls, 1, 2)
! call stop_timer_cpu("cshift_ver3")

call start_timer_cpu("cshift_ver4")
do iptcl=1,nspace
    call progress(iptcl, nspace)
    ptcllim = plim(iptcl)
    do irot=1,nrots
        rotlim = rolim4corr(irot)
        ! calculate hadamard product here
        ! send the below array to gpu
        ! pfts_ptcls(ptcllim(1):ptcllim(2),rotlim(1):rotlim(2),:)
        
    end do
end do
! one last call to get the original order back
pfts_ptcls = cshift(pfts_ptcls, 1, 1)
call stop_timer_cpu("cshift_ver4")
! call start_timer_cpu("cshift_ver5")
! do irot=1,nrots
!     call progress(irot, nrots)
!     lim = plim4corr(irot)
!     do iptcl=1,nspace
!         if( iptcl /= 1 ) pfts_ptcls = cshift(pfts_ptcls, 1, 1)
!         pfts_ptcls_transf(:,:,:) = pfts_ptcls(:,lim(1):lim(2),:)
!         ! calculate hadamard product here
!     end do
! end do
! ! one last call to get the original order back
! pfts_ptcls = cshift(pfts_ptcls, 1, 1)
! call stop_timer_cpu("cshift_ver5")

! mysh(1) = 1
! mysh(2) = 2
! call start_timer_cpu("cshift_ver6")
! do ptclsh=0,nspace-1
!     iptcl = ptclsh+1
!     call progress(iptcl, nspace)
!     do rotsh=0,nrots-1
!         irot = rotsh+1
!
!
!
!         pfts_ptcls = cshift(pfts_ptcls, shift=mysh)
!         ! calculate hadamard product here
!     end do
! end do
! call stop_timer_cpu("cshift_ver6")

call stop_Alltimers_cpu()

contains
    
!     function plim_fake( iptcl ) result( lim )
!         integer, intent(in) :: iptcl
!         integer :: lim(2)
!         lim(1) = (iptcl-fromp)*NROTS_FAKE+1
!         lim(2) = lim(1)+NROTS_FAKE-1
!     end function

!     function plim( iptcl ) result( lim )
!         integer, intent(in) :: iptcl
!         integer :: lim(2)
!         lim(1) = (iptcl-fromp)*nrots+1
!         lim(2) = lim(1)+nrots-1
!     end function
    
    function rolim( r ) result( lim )
        integer, intent(in) :: r
        integer :: lim(2)
        lim(1) = r
        lim(2) = lim(1)+ptclsz-1
    end function
    
    function rolim4corr( r ) result( lim )
        integer, intent(in) :: r
        integer :: lim(2)
        lim(1) = r
        lim(2) = lim(1)+refsz-1
    end function
    
    function plim( iptcl ) result( lim )
        integer, intent(in) :: iptcl
        integer :: lim(2)
        lim(1) = iptcl
        lim(2) = lim(1)+nspace-1
    end function

!     function plim4corr( iptcl ) result( lim )
!         integer, intent(in) :: iptcl
!         integer :: lim(2)
!         lim(1) = (iptcl-fromp)*nrots+1
!         lim(2) = lim(1)+refsz-1
!     end function
    
    subroutine print_fake_pft
        integer :: i, j
        do i=1,NREFS_FAKE
            do j=1,NROTS_FAKE
                write(*,*) fake_pft(i,j,:)
            end do
        end do
        write(*,*) '***********************'
    end subroutine
    
!     subroutine cshift_ptcl_fwd
!         fake_pft = cshift(fake_pft, NROTS_FAKE, 1)
!     end subroutine
!
!     subroutine cshift_ptcl_bwd
!         fake_pft = cshift(fake_pft, -NROTS_FAKE, 1)
!     end subroutine
!
!     subroutine rotate_ptcls_fwd
!         integer :: lim(2), iptcl
!         do iptcl=1,nspace
!             lim = plim(iptcl)
!             pfts_ptcls(lim(1):lim(2),:) = cshift(pfts_ptcls(lim(1):lim(2),:), 1, 1)
!         end do
!     end subroutine
!
!     subroutine rotate_fake_ptcls_fwd
!         integer :: lim(2), iptcl
!         do iptcl=1,NPTCLS
!             lim = plim_fake(iptcl)
!             fake_pft(lim(1):lim(2),:) = cshift(fake_pft(lim(1):lim(2),:), 1, 1)
!         end do
!     end subroutine
!
!     subroutine rotate_fake_ptcls_bwd
!         integer :: lim(2), iptcl
!         do iptcl=1,NPTCLS
!             lim = plim_fake(iptcl)
!             fake_pft(lim(1):lim(2),:) = cshift(fake_pft(lim(1):lim(2),:), -1, 1)
!         end do
!     end subroutine

end program simple_test_cshift