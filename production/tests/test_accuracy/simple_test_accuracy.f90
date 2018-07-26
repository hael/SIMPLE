

program test_expm1
  include 'simple_lib.f08'
  use, intrinsic :: iso_c_binding, only: c_float, c_double
  implicit none
  interface
     function expm1(x) bind(c,name='expm1')
       use, intrinsic :: iso_c_binding, only: c_double
       implicit none
       real(c_double), intent(in) :: x
       real(c_double) :: expm1
     end function expm1

     function expm1f(x) bind(c,name='expm1f')
       use, intrinsic :: iso_c_binding, only: c_float
       implicit none
       real(c_float), intent(in) :: x
       real(c_float) :: expm1f
     end function expm1f
     function sqrtf(x) bind(c, name='sqrtf')
       use, intrinsic :: iso_c_binding, only: c_float
       implicit none
       real(c_float), intent(in) :: x
       real(c_float) :: sqrtf
     end function sqrtf
     function sinhf(x) bind(c, name='sinhf')
       use, intrinsic :: iso_c_binding, only: c_float
       implicit none
       real(c_float), intent(in) :: x
       real(c_float) :: sinhf
     end function sinhf
  end interface

  real(sp), parameter ::  SPECIAL_ENCODING_MASK32 =REAL(Z'60000000',sp)
  real(sp), parameter ::  LARGE_COEFF_MASK32      =REAL(Z'007fffff',sp)
  real(sp), parameter ::  LARGE_COEFF_HIGH_BIT32  =REAL(Z'00800000',sp)
  real(sp), parameter ::  SMALL_COEFF_MASK32      =REAL(Z'001fffff',sp)
  real(sp), parameter ::  EXPONENT_MASK32         =REAL(Z'ff000000',sp)
  real(sp), parameter ::  LARGEST_BID32           =REAL(Z'77f8967f',sp)
  real(sp), parameter ::  SNAN_MASK32             =REAL(Z'7e000000',sp)
  real(sp), parameter ::  INFINITY_MASK32 = REAL(Z'78000000',sp)
  real(sp), parameter ::  NAN_MASK32 = REAL(Z'7c000000',sp)
  integer(dp):: i, t1,n, j
  integer(dp),parameter:: nsize=100000, jsize=5000
  real :: ansme,ansm1c,ansintrinsic
  real(kind=4) :: vector_arg(nsize), sarg(nsize), tmp(nsize)
  real(kind=8) :: tsum(8), dsum(6), dmin(6),dmax(6)

  print *,">>> TEST ACCURACY"
  print *,"  Smallest distance between two sp reals: ", spacing(1.0_sp)
  print *,"  Smallest distance between two dp reals: ", spacing(1.0_dp)
  print *,"  Reciprocal relative distance between two sp reals: ", rrspacing(1.0_sp)
  print *,"  Reciprocal relative distance between two dp reals: ", rrspacing(1.0_dp)
  print *,"  Epsilon for sp reals: ", spacing(1.0_sp)
  print *,"  Epsilon for dp reals: ", spacing(1.0_dp)

  call init_random_seed()
  t1=tic()
  dsum = 0.0
  tsum = 0.0
  dmin = 1.0e20
  dmax = 1.0e-20
  !    do i=1,nsize
  call random_number(vector_arg)
  !     end do
  sarg=vector_arg*5.75+3*TINY
  vector_arg = vector_arg*8.+3.7
  !   print *," Gen args ", toc()
  !   call sleep(1)
  do j =1,jsize
     t1=tic()
     do i=1,nsize
        tmp(i) = exp( vector_arg(i) )/ sqrt(vector_arg(i)) !exp( vector_arg )/ sqrt( vector_arg )
     end do
     tsum(1)=tsum(1)+toc(t1)


     tmp=0.
     t1=tic()
     do i=1,nsize
        tmp(i) = (expm1f( vector_arg(i) ) +1.) / sqrtf( vector_arg(i)) !(vector_arg(i)) )!* InvSqrt( vector_arg(i) )! / sqrtf0( vector_arg(i) )
     end do
     ! tmp = (expm1f( real(vector_arg)) +1. )/ sqrt_flt0( real(vector_arg) )
     tsum(2)=tsum(2)+toc(t1)


     tmp=0.
     t1=tic()
     do i=1,nsize
        tmp(i) =  (exp(vector_arg(i))) * InvSqrt( vector_arg(i)) 
     end do
     !tmp = InvSqrt( real(vector_arg) ) * (expm1_flt_me(real(vector_arg)) +1.)
     tsum(3)=tsum(3)+toc(t1)

     tmp=0.
     t1=tic()
     tmp = (expm1_flt_me((vector_arg)) +1.) / sqrt(vector_arg)
     tsum(7)=tsum(7)+toc(t1)

     tmp=0.
     t1=tic()
     tmp =  (expm1_flt_me(vector_arg) +1.) * InvSqrt( vector_arg ) 
     tsum(8)=tsum(8)+toc(t1)

     tmp=0.
     t1=tic()
     do i=1,nsize
        tmp(i) = sinh( sarg(i) ) / sarg(i)
     end do
     tsum(4)=tsum(4)+toc(t1)
     t1=tic()
     do i=1,nsize
        tmp(i) = sinhf(real(sarg(i))) / sarg(i)
     end do
     !tmp = sinhf(sarg)/sarg
     tsum(5)=tsum(5)+toc(t1)
     t1=tic()
     do i=1,nsize
        tmp(i) = sinch_flt(sarg(i) )
     end do
     !        tmp = sinch_flt(sarg)
     tsum(6)=tsum(6)+toc(t1)
     !        print *, 'arg ',vector_arg(j),'| exp-1 ', exp(vector_arg(j)), '|  my expm1 ',expm1fme(vector_arg(j))+1.
  end do
  print*,"Timing calculations     :  (Time s)            (Speed up)"
  print*,"  exp(x) / sqrt(x)      ", tsum(1)/real(jsize)
  print*,"  expm1(x)/sqrt_flt0(x) ", tsum(2)/real(jsize),"   ", (tsum(1)-tsum(2))/tsum(1) ,  tsum(1)/tsum(2)
  print*,"  expm1f(x)*InvSqrt(x)  ", tsum(3)/real(jsize),"   ", (tsum(1)-tsum(3))/tsum(1) , tsum(1)/tsum(3)
  print*,"  expm1_flt_me(x)*InvSqrt(x)  ", tsum(8)/real(jsize),"   ", (tsum(1)-tsum(8))/tsum(1) , tsum(1)/tsum(8)
  print*,"  expm1_flt_me(x)/sqrt(x)  ", tsum(7)/real(jsize),"   ", (tsum(1)-tsum(7))/tsum(1) , tsum(1)/tsum(7)

  print*,''
  print*,''
  print*,"  sinh(x)/x              ", tsum(4)/real(jsize)
  print*,"  sinhf(x)/x             ", tsum(5)/real(jsize), "   ", (tsum(4)-tsum(5))/tsum(4), tsum(4)/tsum(5)
  print*,"  sinch(x)               ", tsum(6)/real(jsize), "   ", (tsum(4)-tsum(6))/tsum(4) , tsum(4)/tsum(5)
  tsum=0.
  dsum=0.
  print*,''
  print*,''
  print*,''

  !call sleep(1)
  do i=1,nsize/100

     ansintrinsic = exp( vector_arg(i) ) / sqrt( vector_arg(i) ) !!  exp(vector_arg(i) )
     ansm1c  = exp( vector_arg(i) )  / sqrt_flt0( real(vector_arg(i)) )
     ansme   =  InvSqrt( real( vector_arg(i) ) ) *  exp( vector_arg(i) )
     dsum(1) = dsum(1)+ (abs(ansintrinsic-ansme))
     if (abs(ansintrinsic-ansme)< dmin(1) )dmin(1)=abs(ansintrinsic-ansme)
     if (abs(ansintrinsic-ansme)> dmax(1) )dmax(1)=abs(ansintrinsic-ansme)
     dsum(2) = dsum(2)+ (abs(ansintrinsic-ansm1c))
     if (abs(ansintrinsic-ansm1c)< dmin(2) )dmin(2)=abs(ansintrinsic-ansm1c)
     if (abs(ansintrinsic-ansm1c)> dmax(2) )dmax(2)=abs(ansintrinsic-ansm1c)
     ansintrinsic  = sinh( sarg(i) ) / sarg(i)
     ansme   = sinch_dbl(REAL(sarg(i),8))
     dsum(3) = dsum(3)+ (abs(ansintrinsic-ansme))
     if (abs(ansintrinsic-ansme)< dmin(3) )dmin(3)=abs(ansintrinsic-ansme)
     if (abs(ansintrinsic-ansme)> dmax(3) )dmax(3)=abs(ansintrinsic-ansme)
     ansm1c  = sinh_flt_me(real(sarg(i)) )/ sarg(i)
     dsum(4) = dsum(4)+ (abs(ansintrinsic-ansm1c))
     if (abs(ansintrinsic-ansm1c)< dmin(4) )dmin(4)=abs(ansintrinsic-ansm1c)
     if (abs(ansintrinsic-ansm1c)> dmax(4) )dmax(4)=abs(ansintrinsic-ansm1c)
     ! if(.not.tolerence(ansintrinsic,ansme,1.e-6)) then
     !     print *,' Error in expm1 ', abs(ansintrinsic-ansme)
     !     stop
     ! end if

  end do
  print*,"Error calculation      (Mean)                 (Min)             (Max)"
  print*,"  expm1f(x)*InvSqrt(x)  ", dsum(1)/REAL(nsize),  dmin(1),dmax(1)
  print*,"  exp*InvSqrt           ", dsum(2)/REAL(nsize),  dmin(2),dmax(2)

  print*,"  sinhf(x)/x            ", dsum(4)/REAL(nsize),  dmin(4),dmax(4)
  print*,"  sinch(x)              ", dsum(3)/REAL(nsize),  dmin(3),dmax(3)


contains
  subroutine init_random_seed()

    integer, allocatable :: seed(:)
    integer              :: i, n, un, istat, dt(8), pid, t(2), s
    integer(8)           :: count, tms

    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(count)
       if (count /= 0) then
          t = transfer(count, t)
       else
          call date_and_time(values=dt)
          tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24 * 60 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
          t = transfer(tms, t)
       end if
       s = ieor(t(1), t(2))
       pid = INT(getpid()) + 1099279 ! Add a prime
       s = ieor(s, pid)
       if (n >= 3) then
          seed(1) = t(1) + 36269
          seed(2) = t(2) + 72551
          seed(3) = pid
          if (n > 3) then
             seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
          end if
       else
          seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
       end if
    end if
    call random_seed(put=seed)
  end subroutine init_random_seed

  ! logical pure function is_nan(x)
  !     real(sp), intent(in) :: x

  !     is_nan=.false.
  !     if((x.and.NAN_MASK32) == NAN_MASK32) is_nan=.true.
  ! end function is_nan
  ! logical pure function is_inf(x)
  !     real(sp), intent(in) :: x


  !     is_inf=.false.
  !     if( and(x , INFINITY_MASK32) == INFINITY_MASK32) is_inf=.true.
  ! end function is_inf

  elemental function expm1_flt_me(x) result(y)
    real(sp), intent(in) :: x
    real(sp)             :: y
    integer(dp) :: k
    real(dp) :: absx,hxs,hfs,hfx,hi,r1,lo,t,c,e
    real(sp), parameter :: Othreshold  = DBLE(7.09782712893383973096e+02)  !!DBLE( Z'40862E42FEFA39EF') !!
    real(sp), parameter :: Ln2X56      = DBLE(3.88162421113569373274e+01)  !!DBLE( Z'4043687a9f1af2b1') !!
    real(sp), parameter :: Ln2HalfX3   = DBLE(1.03972077083991796413e+00)  !!DBLE( Z'3ff0a2b23f3bab73') !!
    real(sp), parameter :: Ln2Half     = DBLE(3.46573590279972654709e-01)  !!DBLE( Z'3fd62e42fefa39ef') !!
    real(sp), parameter :: Ln2Hi       = DBLE(6.93147180369123816490e-01)  !!DBLE( Z'3fe62e42fee00000') !!
    real(sp), parameter :: Ln2Lo       = DBLE(1.90821492927058770002e-10)  !!DBLE( Z'3dea39ef35793c76') !!
    real(sp), parameter :: InvLn2      = DBLE(1.44269504088896338700e+00)  !!DBLE( Z'3ff71547652b82fe') !!
    ! real(sp), parameter :: TinyDble    = DBLE( Z'3c90000000000000')  !! 2**-54 = Z'3c90000000000000' !!1.0 / (1 << 54)
    ! scaled coefficients related to expm1
    real(sp), parameter :: Q1 = DBLE(Z'BFA11111111110F4')  !!  -3.33333333333331316428e-02)  !!
    real(sp), parameter :: Q2 = DBLE(Z'3F5A01A019FE5585')  !!  1.58730158725481460165e-03)   !!
    real(sp), parameter :: Q3 = DBLE(Z'BF14CE199EAADBB7')  !!  -7.93650757867487942473e-05)  !!
    real(sp), parameter :: Q4 = DBLE(Z'3ED0CFCA86E65239')  !!  4.00821782732936239552e-06)   !!
    real(sp), parameter :: Q5 = DBLE(Z'BE8AFDB76E09C32D')  !!  -2.01099218183624371326e-07)  !!
    real(sp), parameter :: dxf = 5e-4
    real(dp), parameter :: dxd = 2e-8
    logical :: sign

    !! special cases
    !  if( x/=x )then
    !         y=x
    !     return
    ! else if ( x > HUGE(Q1)) then
    !      y=-1.
    !      return
    !  end if

    ! Case 0: x == 0
    ! y=1.
    ! return
    ! Case 1: |x|<dx,
    ! y=1.+x/2.
    ! return  Choose dx empirically for your numerical precision and data type. For double, it should be about 2e-8. For float, it's about 5e-4.
    ! Case else: return
    ! expm1(x)/x.

    absx = x
    sign = .false.
    ! if ( x < 0 )then
    !     absx = -absx
    !     sign = .true.
    ! end if

    ! !! filter out huge argument
    ! if (absx >= Ln2X56) then !! if |x| >= 56 * ln2
    !     if ( sign )then
    !         y=-1.0
    !         return  !! x < -56*ln2, return -1
    !     end if
    !     if (absx >= Othreshold) then !! if |x| >= 709.78...
    !         y=REAL(INFINITY_MASK32)
    !         return
    !     end if
    ! end if

    !! argument reduction

    if (absx > Ln2Half) then !! if  |x| > 0.5 * ln2

       if (absx < Ln2HalfX3) then !! and |x| < 1.5 * ln2
          !   if (.not.sign)then
          hi = x - Ln2Hi
          lo = Ln2Lo
          k = 1
          ! else
          !     hi = x + Ln2Hi
          !     lo = -Ln2Lo
          !     k = -1
          ! end if
       else
          ! if (.not.sign)then !
          k = int(InvLn2*x + 0.5)
          ! else
          !     k = int(InvLn2*x - 0.5)
          ! end if
          t = k
          hi = x - t*Ln2Hi !! t * Ln2Hi is exact here
          lo = t * Ln2Lo
       end if
       c = hi - lo
       c = (hi - c) - lo
    else if (absx < dxd) then !! when |x| < 2**-54, return x
       y=x
       return
    else
       k = 0
    end if

    !! x is now in primary range
    hfx = 0.5 * x
    hxs = x * hfx
    r1 = 1 + hxs*(Q1+hxs*(Q2+hxs*(Q3+hxs*(Q4+hxs*Q5))))
    t = 3 - r1*hfx
    e = hxs * ((r1 - t) / (6.0 - x*t))
    if( k /= 0 )then
       e = (x*(e-c) - c)
       e = e-hxs
       if(k== -1)then
          y=0.5*(x-e) - 0.5
          return
       else if(k==  1)then
          if (x < -0.25) then
             y = -2 * (e - (x + 0.5))
          else
             y=1 + 2*(x-e)
          end if
          return
       else if( k <= -2 .or. k > 56)then !! suffice to return exp(x)-1
          y = 1 - (e - x)
          y = y+2**k -1 !!Float64frombits(Float64bits(y) + uint64(k)<<52) !! add k to y's exponent
          return
       end if
       if (k < 20) then
          t = 1-2**(-k)  !!Float64frombits(Z'3ff0000000000000' - (Z'20000000000000' >> uint(k)))
          y = t - (e - x)
          y = y+ 2**k !! Float64frombits(Float64bits(y) + uint64(k)<<52) !! add k to y's exponent
          return
       end if
       t = 2**(-k) !! Float64frombits(uint64(Z'3ff-k) << 52)
       y = x - (e + t)

       y = 1+ y+ 2**k !! Float64frombits(Float64bits(y) + uint64(k)<<52) !! add k to y's exponent
       return
    end if
    y=x - (x*e - hxs) !! c is 0

  end function expm1_flt_me

  elemental function expm1_dbl_me(x) result(y)
    real(sp), intent(in) :: x
    real(sp)             :: y
    integer(dp) :: k
    real(dp) :: absx,hxs,hfs,hfx,hi,r1,lo,t,c,e
    real(sp), parameter :: Othreshold  = DBLE(7.09782712893383973096e+02)  !!DBLE( Z'40862E42FEFA39EF') !!
    real(sp), parameter :: Ln2X56      = DBLE(3.88162421113569373274e+01)  !!DBLE( Z'4043687a9f1af2b1') !!
    real(sp), parameter :: Ln2HalfX3   = DBLE(1.03972077083991796413e+00)  !!DBLE( Z'3ff0a2b23f3bab73') !!
    real(sp), parameter :: Ln2Half     = DBLE(3.46573590279972654709e-01)  !!DBLE( Z'3fd62e42fefa39ef') !!
    real(sp), parameter :: Ln2Hi       = DBLE(6.93147180369123816490e-01)  !!DBLE( Z'3fe62e42fee00000') !!
    real(sp), parameter :: Ln2Lo       = DBLE(1.90821492927058770002e-10)  !!DBLE( Z'3dea39ef35793c76') !!
    real(sp), parameter :: InvLn2      = DBLE(1.44269504088896338700e+00)  !!DBLE( Z'3ff71547652b82fe') !!
    ! real(sp), parameter :: TinyDble    = DBLE( Z'3c90000000000000')  !! 2**-54 = Z'3c90000000000000' !!1.0 / (1 << 54)
    ! scaled coefficients related to expm1
    real(sp), parameter :: Q1 = DBLE(Z'BFA11111111110F4')  !!  -3.33333333333331316428e-02)  !!
    real(sp), parameter :: Q2 = DBLE(Z'3F5A01A019FE5585')  !!  1.58730158725481460165e-03)   !!
    real(sp), parameter :: Q3 = DBLE(Z'BF14CE199EAADBB7')  !!  -7.93650757867487942473e-05)  !!
    real(sp), parameter :: Q4 = DBLE(Z'3ED0CFCA86E65239')  !!  4.00821782732936239552e-06)   !!
    real(sp), parameter :: Q5 = DBLE(Z'BE8AFDB76E09C32D')  !!  -2.01099218183624371326e-07)  !!
    real(sp), parameter :: dxf = 5e-4
    real(dp), parameter :: dxd = 2e-8
    logical :: sign

    !! special cases
    !  if( x/=x )then
    !         y=x
    !     return
    ! else if ( x > HUGE(Q1)) then
    !      y=-1.
    !      return
    !  end if

    ! Case 0: x == 0
    ! y=1.
    ! return
    ! Case 1: |x|<dx,
    ! y=1.+x/2.
    ! return  Choose δδ empirically for your numerical precision and data type. For double, it should be about 2e-8. For float, it's about 5e-4.
    ! Case else: return
    ! expm1(x)/x.

    absx = x
    sign = .false.
    ! if ( x < 0 )then
    !     absx = -absx
    !     sign = .true.
    ! end if

    ! !! filter out huge argument
    ! if (absx >= Ln2X56) then !! if |x| >= 56 * ln2
    !     if ( sign )then
    !         y=-1.0
    !         return  !! x < -56*ln2, return -1
    !     end if
    !     if (absx >= Othreshold) then !! if |x| >= 709.78...
    !         y=REAL(INFINITY_MASK32)
    !         return
    !     end if
    ! end if

    !! argument reduction

    if (absx > Ln2Half) then !! if  |x| > 0.5 * ln2

       if (absx < Ln2HalfX3) then !! and |x| < 1.5 * ln2
          !   if (.not.sign)then
          hi = x - Ln2Hi
          lo = Ln2Lo
          k = 1
          ! else
          !     hi = x + Ln2Hi
          !     lo = -Ln2Lo
          !     k = -1
          ! end if
       else
          ! if (.not.sign)then !
          k = int(InvLn2*x + 0.5)
          ! else
          !     k = int(InvLn2*x - 0.5)
          ! end if
          t = k
          hi = x - t*Ln2Hi !! t * Ln2Hi is exact here
          lo = t * Ln2Lo
       end if
       c = hi - lo
       c = (hi - c) - lo
    else if (absx < dxd) then !! when |x| < 2**-54, return x
       y=x
       return
    else
       k = 0
    end if

    !! x is now in primary range
    hfx = 0.5 * x
    hxs = x * hfx
    r1 = 1 + hxs*(Q1+hxs*(Q2+hxs*(Q3+hxs*(Q4+hxs*Q5))))
    t = 3 - r1*hfx
    e = hxs * ((r1 - t) / (6.0 - x*t))
    if( k /= 0 )then
       e = (x*(e-c) - c)
       e = e-hxs
       if(k== -1)then
          y=0.5*(x-e) - 0.5
          return
       else if(k==  1)then
          if (x < -0.25) then
             y = -2 * (e - (x + 0.5))
          else
             y=1 + 2*(x-e)
          end if
          return
       else if( k <= -2 .or. k > 56)then !! suffice to return exp(x)-1
          y = 1 - (e - x)
          y = y+2**k -1 !!Float64frombits(Float64bits(y) + uint64(k)<<52) !! add k to y's exponent
          return
       end if
       if (k < 20) then
          t = 1-2**(-k)  !!Float64frombits(Z'3ff0000000000000' - (Z'20000000000000' >> uint(k)))
          y = t - (e - x)
          y = y+ 2**k !! Float64frombits(Float64bits(y) + uint64(k)<<52) !! add k to y's exponent
          return
       end if
       t = 2**(-k) !! Float64frombits(uint64(Z'3ff-k) << 52)
       y = x - (e + t)

       y = 1+ y+ 2**k !! Float64frombits(Float64bits(y) + uint64(k)<<52) !! add k to y's exponent
       return
    end if
    y=x - (x*e - hxs) !! c is 0

  end function expm1_dbl_me



  logical function tolerence (a1, b1, e1 )
    real, intent(in) :: a1,b1,e1
    real :: d, esp

    !! Multiplying by e here can underflow denormal values to zero.
    !! Check a==b so that at least if a and b are small and identical
    !! we say they match.
    if (a1 == b1) then
       tolerence = .true.
       return
    end if
    d = a1 - b1
    if (d < 0) then
       d = -d
    end if

    !! note: b is correct (expected) value, a is actual value.
    !! make error tolerance a fraction of b, not a.
    if (b1 /= 0.0) then

       esp = e1 * b1
       if ( esp < 0.0 ) then
          esp = -esp

       end if
       tolerence= (d < esp)
    end if
  end function tolerence

  elemental pure function sinch_dbl(xin) result (y)
    !! The coefficients are #2029 from Hart & Cheney. (20.36D)
    real(dp), intent(in) :: xin
    real(dp), parameter  :: P0 = -0.6307673640497716991184787251d+6,&
         P1 = -0.8991272022039509355398013511d+05, &
         P2 = -0.2894211355989563807284660366d+04, &
         P3 = -0.2630563213397497062819489000d+02, &
         Q0 = -0.6307673640497716991212077277d+06, &
         Q1 =  0.1521517378790019070696485176d+05, &
         Q2 = -0.1736789535582336995334509110d+03
    !logical ::sign
    real(sp) :: y,x,xsq

    x = dble(xin)
    !sign = .false.
    ! if (x < 0. ) then
    !     x = -x
    !     sign = .true.
    ! end if
    ! if (x > 21.) then
    !     y = exp(x) / (x*2.)
    ! else
    if (x > 0.5) then
       y = (exp(x) - exp(-x)) / (x*2.)
    else
       xsq = x * x
       y = (((P3*xsq+P2)*xsq+P1)*xsq + P0)
       y = y / (((xsq+Q2)*xsq+Q1)*xsq + Q0)
    end if

    ! if (sign) then
    !     y = -y
    ! end if

  end function sinch_dbl
  elemental pure function sinch_flt(x) result (y)
    !! The coefficients are #2029 from Hart & Cheney. (20.36D)
    real(sp), intent(in) :: x
    real(sp), parameter  :: P0 = -0.6307673640497716991184787251d+6,&
         P1 = -0.8991272022039509355398013511d+05, &
         P2 = -0.2894211355989563807284660366d+04, &
         P3 = -0.2630563213397497062819489000d+02, &
         Q0 = -0.6307673640497716991212077277d+06, &
         Q1 =  0.1521517378790019070696485176d+05, &
         Q2 = -0.1736789535582336995334509110d+03
    !logical ::sign
    real(sp) :: y,xsq


    !sign = .false.
    ! if (x < 0. ) then
    !     x = -x
    !     sign = .true.
    ! end if
    ! if (x > 21.) then
    !     y = exp(x) / (x*2.)
    ! else
    if (x > 0.5) then
       y = (exp(x) - exp(-x)) / (x*2.)
    else
       xsq = x * x
       y = (((P3*xsq+P2)*xsq+P1)*xsq + P0)
       y = y / (((xsq+Q2)*xsq+Q1)*xsq + Q0)
    end if

    ! if (sign) then
    !     y = -y
    ! end if

  end function sinch_flt


  elemental pure function sinh_flt_me(xin) result (y)
    !! The coefficients are #2029 from Hart & Cheney. (20.36D)
    real(sp), intent(in) :: xin
    real(dp), parameter  :: P0 = -0.6307673640497716991184787251d+6,&
         P1 = -0.8991272022039509355398013511d+05, &
         P2 = -0.2894211355989563807284660366d+04, &
         P3 = -0.2630563213397497062819489000d+02, &
         Q0 = -0.6307673640497716991212077277d+06, &
         Q1 =  0.1521517378790019070696485176d+05, &
         Q2 = -0.1736789535582336995334509110d+03
    !logical ::sign
    real(dp) :: y,x,xsq
    x=xin
    ! sign = .false.
    ! if (x < 0. ) then
    !     x = -x
    !     sign = .true.
    ! end if
    ! if (x > 21.) then
    !     y = exp(x) / 2.
    ! else
    if (x > 0.5) then
       y = (exp(x) - exp(-x)) / 2.
    else
       xsq = x * x
       y = (((P3*xsq+P2)*xsq+P1)*xsq + P0) * x
       y = y / (((xsq+Q2)*xsq+Q1)*xsq + Q0)
    end if

    ! if (sign) then
    !     y = -y
    ! end if

  end function sinh_flt_me

  logical function isinf(x)
    ! * (is Infinity?)
    ! * Return .TRUE. if x is Infinity, and
    ! * .FALSE. otherwise, without causing a
    ! * floating-point trap.
    ! * (09-Mar-1994)
    real ::x,wf
    integer ::wi
    equivalence (wf, wi)
    wf = x
    ! Inf has maximum exponent, and zero fraction
    isinf = (and(rshift(wi,23),255) .eq. 255) .and. (and(wi,8388607) .eq. 0)

  end function isinf
  logical function isnan(x)
    ! * (is NaN?)
    ! * Return .TRUE. if x is a NaN, and
    ! * .FALSE. otherwise, without causing a
    ! * floating-point trap.
    ! * (09-Mar-1994)
    real ::x,wf
    integer ::wi
    equivalence (wf, wi)
    wf = x
    ! NaN has maximum exponent, and non-zero fraction
    isnan = (and(rshift(wi,23),255) .eq. 255) .and. (and(wi,8388607) .ne. 0)
  end function isnan
  elemental pure logical function isden(x)
    ! * (is denormalized?)
    ! * Return .TRUE. if x is denormalized, and
    ! * .FALSE. otherwise, without causing a
    ! * floating-point trap.
    ! * (09-Mar-1994)
    real(4), intent(in) :: x
    integer(4) :: wi
    wi = transfer(x,wi)
    ! * Denorm has minimum exponent, and non-zero
    ! * fraction
    isden = (iand(shiftr(wi,23),255) .eq. 0) .and. (iand(wi,8388607) .ne. 0)
  end function isden
  elemental pure real function setxp (x,n)
    ! (Set exponent of x)
    ! [14-Nov-1990]
    real , intent(in) :: x
    integer , intent(in):: n
    integer :: wi
    wi=transfer(x,wi)

    ! Zero exponent field
    wi = iand(wi,(z'807fffff'))

    ! Or in the new exponent field
    wi = ior(wi,shiftl(iand(126 + n,(z'ff')),23))

    setxp = transfer(wi,1.0)

  end function setxp

  elemental pure integer function intxp(x)
    ! (Return unbiased exponent of x)
    ! [14-Nov-1990]

    real , intent(in) :: x
    integer :: wi

    ! Force storage overlay so we can twiddle bits
    wi=transfer(x,wi)


    ! Extract the exponent field and unbias
    intxp = iand(shiftr(wi,23),INT(z'ff')) - 126

  end function intxp
  elemental pure real function adx(x,n)
    ! (Adjust exponent of x)
    ! [14-Nov-1990]
    real , intent(in) :: x
    integer , intent(in):: n
    integer :: wi, olde
    wi=transfer(x,wi)
    ! Extract old exponent
    olde = iand(shiftr(wi,23),INT(z'ff'))

    ! Increment old exponent
    olde = iand(olde + n, INT(z'ff'))

    ! Zero the exponent fieldd
    wi = iand(wi,(z'807fffff'))

    ! Or in the new exponent
    wi = ior(wi,shiftl(olde,23))

    adx = transfer(wi,1.0)

  end function adx

  elemental pure real function sqrt_flt0(x)
    !  Cody-Waite implementation of sqrt(x)
    real(sp), intent(in):: x
    real(sp) :: y,z,f,xx
    integer  :: exponent, nbias
    real, parameter :: denmax = REAL(z'007fffff') ! = maximum positive denormalized number
    real, parameter :: onemme = REAL(z'3f7fffff')! = 1.0 - machine epsilon
    !        real, parameter ::  Inf = 0./0.! REAL(z'7f800000')

    ! isden(xx) = (xx .le. denmax)
    ! isinf(xx) = (xx .ge. Inf)
    ! isnan(xx) = (xx .ne. xx)
    if (isden(x)) then
       ! scale x by 2**24 (this is exact)
       xx = x * 16777216.0
       nbias = -24
    else
       xx = x
       nbias = 0
    end if
    exponent = intxp(xx)
    f = setxp(xx, 0)
    y = 0.41731 + 0.59016 *f
    !    y1 = 0.5 * (y0 + f /y0)

    !    y2 = 0.5 * (y1 + f /y1)

    z = (y + f /y)
    y = 0.5 * z
    y = 0.25 * z + f /z
    if (MOD(exponent,2).ne.0)then
       y = y * 0.707106781186547524400844362104E+00
       y = max(y,0.5)
       exponent=exponent+1
    end if
    y = min(y,onemme)

    ! Insert exponent to undo range reduction.
    sqrt_flt0 = setxp(y,(exponent + nbias)/2)

  end function sqrt_flt0

  elemental pure function InvSqrt( xin ) result(y)
    real, intent(in) :: xin
    real(4) :: x,xhalf,y
    integer(4) :: i
    integer(4), parameter :: QUAKEMAGIC = int(Z'5F375A86',4)!! INT(Z'5FE6EB50C7B537A9',4)
    !!int(Z'5F375A86',4) !!INT(Z'5f3759df') !! 1597463007_sp 1597463007_sp5F375A86,
    !!0x5FE6EB50C7B537A9
    equivalence(x,i)            !! store floating-point bits in integer
    x=xin
    xhalf = 0.5 * x
    i = QUAKEMAGIC - SHIFTR(i,1) !! initial guess for Newton's method
    !x = transfer(i,x)           !! convert new bits into float
    y  = x * (1.5 - ( xhalf * x * x ) )     !! One round of Newton's method
    y  = y * (1.5 - ( xhalf * y * y ) )
  end function InvSqrt
end program test_expm1
