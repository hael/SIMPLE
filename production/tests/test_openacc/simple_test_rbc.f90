!! Dynamic equilibrium model http://www.kellogg.northwestern.edu/rc/workshops/mp2.html
!! 
module simple_rbc
  use omp_lib     ! For timing
#ifndef __INTEL_COMPILER
  use openacc
#endif
  implicit none
  private

  real, parameter :: beta = 0.984 ! discount rate
  real, parameter :: eta = 2.!Risk aversion parameter
  real, parameter ::alpha = 0.35!Technology parameter
  real, parameter :: delta = 0.01!!Depreciation rate
  real, parameter :: rho = 0.95!Tech. shock persistence
  real, parameter :: sigma = 0.005!Tech. shock st. dev.
  real, parameter :: zmin=-0.0480384!Grid for productivity z
  real, parameter :: zmax=0.0480384!
 real, parameter :: tol = 1e-4
  integer, parameter :: nz = 4! Grid size
  integer, parameter :: nk=4800!! Number of points for capital
  real :: t_tran_z(nz,nz), tran_z(nz,nz)
  real :: zgrid(nz), kgrid(nk)
  real :: kmax, kmin, kstep, dif

  real, allocatable :: v(:,:), v0(:,:), ev(:,:), c0(:,:), c(:), r(:), w(:)
  integer :: i, iz, ik, cnt
  logical :: ind(nk)
  real(kind=8) :: start, finish   ! For timing
  real :: tmpmax, c1, zstep

  public:: run_rbc_oacc, run_rbc_omp, run_rbc_serial
contains

  elemental real function uf(c)
    real , intent(in):: c
    ! Define the utility function
    uf = c**(1-eta) / (1-eta)
  end function uf

  subroutine rbc_init
    integer :: iostat
    start = omp_get_wtime()


    allocate( v(nk,nz), v0(nk,nz), ev(nk,nz), c0(nk,nz), stat=iostat)
    if(iostat/=0) stop 'store'

    allocate( c(nk), r(nk), w(nk), stat=iostat)
    if(iostat/=0) stop 'temps'
    v=0.; v0=0.; ev=0.; c0=0.
    zstep = ((zmax-zmin)/(nz-1))
    ! [1 x 4] grid of values for z
    do iz = 1,nz
       zgrid(iz) = exp(zmin + (iz-1)*zstep )
    end do
    ! [4 x 4] Markov transition matrix of z
    tran_z(:,1) = (/  0.996757 , 0.00324265 , 0.         , 0. /)
    tran_z(:,2) = (/      0.000385933     , 0.998441   , 0.00117336 , 0. /)
    tran_z(:,3) = (/        0.              , 0.00117336 , 0.998441   , 0.000385933 /)
    tran_z(:,4) = (/        0.              , 0.         , 0.00324265 , 0.996757 /)
    t_tran_z = transpose(tran_z)
    kmin = 0.95*(1/(alpha*zgrid(1)))*((1/beta)-1+delta)**(1/(alpha-1));
    kmax = 1.05*(1/(alpha*zgrid(nz)))*((1/beta)-1+delta)**(1/(alpha-1));
    kstep = (kmax - kmin)/(nk-1);
    ! [1 x 4800] grid of possible values of k
    do ik=1,nk
       kgrid(ik) = kmin + (ik-1)* kstep
    end do

    ! Compute initial wealth c0(k,z)
    do iz = 1,nz
       c0(:,iz) = zgrid(iz)*(kgrid**alpha) + (1-delta)*kgrid
    end do

    dif=huge(0.)
    cnt=0

    finish = omp_get_wtime()
    print *,' rbc_init ', finish -start


  end subroutine rbc_init

  subroutine run_rbc_serial
    call rbc_init
    start = omp_get_wtime()
    do while(dif>tol)
       do iz=1,nz
          do ik = 1,nk
             c = c0(ik,iz)-kgrid
             ind = c>0;
             v(ik,iz) = maxval(pack(c,ind)**(1-eta)/(1-eta)+pack(ev(:,iz), ind));
          end do
       end do
       ev = beta*matmul(v,t_tran_z)
       dif = maxval(abs(v-v0))
       v0 = v
       if(mod(cnt,100)==0) write(*,*) cnt, ':', dif
       cnt = cnt+1
    end do

    finish = omp_get_wtime()
    print *,' run_rbc_serial ', finish -start, v(nk,nz)
    deallocate( v, v0, ev, c0, c,r,w)

  end subroutine run_rbc_serial


  subroutine run_rbc_serial1
    call rbc_init
    start = omp_get_wtime()
    do while(dif>tol)

       do iz=1,nz
          do ik = 1,nk
             c = c0(ik,iz)-kgrid
             ind = c>0;
             v(ik,iz) = maxval(pack(uf(c),ind)+pack(ev(:,iz), ind));
          end do
       end do
       ev = beta*matmul(v,t_tran_z)
       dif = maxval(abs(v-v0))
       v0 = v
       if(mod(cnt,100)==0) write(*,*) cnt, ':', dif
       cnt = cnt+1
    end do

    finish = omp_get_wtime()
    print *,' run_rbc_serial ', finish -start, v(nk,nz)

    deallocate( v, v0, ev, c0, c,r,w)

  end subroutine run_rbc_serial1

  subroutine run_rbc_omp
    call rbc_init

    start = omp_get_wtime()
    do while(dif>tol)
       !$omp parallel do collapse(2) default(shared) private(ik,iz,i,tmpmax,c1)
       do ik=1,nk
          do iz = 1,nz
             tmpmax = -huge(0.)
             do i = 1,nk
                c1 = c0(ik,iz) - kgrid(i)
                if(c1<0) exit
                c1 = c1**(1-eta)/(1-eta)+ev(i,iz)
                if(tmpmax<c1) tmpmax = c1
             end do
             v(ik,iz) = tmpmax
          end do
       end do
       !$omp end parallel do
       ev = beta*matmul(v,t_tran_z)
       dif = maxval(abs(v-v0))
       v0 = v
       if(mod(cnt,100)==0) write(*,*) cnt, ':', dif
       cnt = cnt+1
    end do

    finish = omp_get_wtime()
    print *,' run_rbc_omp cnt ', cnt, ' time ', finish -start, 'val ', v(nk,nz)

    deallocate( v, v0, ev, c0, c,r,w)
  end subroutine run_rbc_omp

  subroutine run_rbc_oacc
    call rbc_init

    start = omp_get_wtime()
    do while(dif>tol)
       !$acc parallel
       !$acc loop collapse(2) private(ik,iz,i,tmpmax,c1)
       do ik = 1,nk
          do iz=1,nz;
             tmpmax = -huge(0.)
             do i = 1,nk
                c1 = c0(ik,iz) - kgrid(i)
                if(c1<0) exit
                c1 = c1**(1-eta)/(1-eta)+ev(i,iz)
                if(tmpmax<c1) tmpmax = c1
             end do
             v(ik,iz) = tmpmax
          end do
       end do
       !$acc end loop
       !$acc end parallel
       ev = beta*matmul(v,t_tran_z)
       dif = maxval(abs(v-v0))
       v0 = v
        if(mod(cnt,100)==0) write(*,*) cnt, ':', dif
       cnt = cnt+1
    end do

    finish = omp_get_wtime()
    print *,' run_rbc_oacc ', finish -start, v(nk,nz)
    deallocate( v, v0, ev, c0, c,r,w)
  end subroutine run_rbc_oacc


end module simple_rbc
