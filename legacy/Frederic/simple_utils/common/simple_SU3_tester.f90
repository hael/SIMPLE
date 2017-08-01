!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 9th of October 2013.
!
! Name:
! tester - Various utilities and matrix getter for other modules.
!
! Description:
! tester module provides test code.
! Subroutine su3randomlinks_mul_timing_cpu(elps_tm,t1s,t2s)
!            : test the timers for the CPU 
!*******************************************************************************
!
module simple_SU3_tester

  use simple_defs
  use simple_timing
  use simple_random
  use matrixGetter

  implicit none

contains

  subroutine su3randomlinks_mul_timing_cpu(elps_tm,t1s,t2s)
    implicit none

    !global variables

    double precision,dimension(3)                           :: elps_tm
    double precision,dimension(3)                           :: t1s,t2s

    !local variables

    double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: u1r_rnd,u1i_rnd
    double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: u2r_rnd,u2i_rnd
    double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: u3r,u3i

    !timing variables
    double precision                :: elps_tm_mm
    double precision,dimension(2)   :: s_tm_mm, e_tm_mm
    double precision :: t1_mm,t2_mm

    double precision                :: elps_tm_u1, elps_tm_u2
    double precision,dimension(2)   :: s_tm_u1, e_tm_u1
    double precision,dimension(2)   :: s_tm_u2, e_tm_u2
    double precision :: t1_u1,t2_u1
    double precision :: t1_u2,t2_u2

    integer                         :: ic,jc,kc  !loop indices
    !start of the excution commands

    !initialising the fixed seed for the random number generator
    call init_fixed_random_seed(1234)

    !serial multiplication

    call start_timer_cpu("su3random")
    call cpu_time(t1_u1)
#if defined (LINUX)
    call gettimeofday_c(s_tm_u1)
#endif
    call su3random(u1r_rnd,u1i_rnd)
#if defined (LINUX)
    call gettimeofday_c(e_tm_u1)
    call elapsed_time_c(s_tm_u1,e_tm_u1,elps_tm_u1)
#endif
    call cpu_time(t2_u1)
    call stop_timer_cpu("su3random")

    call cpu_time(t1_u2)
#if defined (LINUX)
    call gettimeofday_c(s_tm_u2)
#endif
    call su3random(u2r_rnd,u2i_rnd)
#if defined (LINUX)
    call gettimeofday_c(e_tm_u2)
    call elapsed_time_c(s_tm_u2,e_tm_u2,elps_tm_u2)
#endif
    call cpu_time(t2_u2)

    call cpu_time(t1_mm)
#if defined (LINUX)
    call gettimeofday_c(s_tm_mm)
#endif
    u3r = 0.0d0
    u3i = 0.0d0
    do ic=1,nc
       do jc=1,nc
          do kc=1,nc
             u3r(:,:,:,:,:,ic,jc) = u3r(:,:,:,:,:,ic,jc) + &
                  ( u1r_rnd(:,:,:,:,:,ic,kc) * u2r_rnd(:,:,:,:,:,kc,jc) - &
                  u1i_rnd(:,:,:,:,:,ic,kc) * u2i_rnd(:,:,:,:,:,kc,jc) )
             u3i(:,:,:,:,:,ic,jc) = u3i(:,:,:,:,:,ic,jc) + &
                  ( u1r_rnd(:,:,:,:,:,ic,kc) * u2i_rnd(:,:,:,:,:,kc,jc) + &
                  u1i_rnd(:,:,:,:,:,ic,kc) * u2r_rnd(:,:,:,:,:,kc,jc) )
          end do
       end do
    end do

#if defined (LINUX)
    call gettimeofday_c(e_tm_mm)
    call elapsed_time_c(s_tm_mm,e_tm_mm,elps_tm_mm)
#endif
    call cpu_time(t2_mm)

    elps_tm(1) = elps_tm_u1
    elps_tm(2) = elps_tm_u2
    elps_tm(3) = elps_tm_mm

    t1s(1) = t1_u1
    t1s(2) = t1_u2
    t1s(3) = t1_mm

    t2s(1) = t2_u1
    t2s(2) = t2_u2
    t2s(3) = t2_mm

    return
  end subroutine su3randomlinks_mul_timing_cpu

  subroutine su3randomlinks_mul_timing(elps_tm,t1s,t2s)
    implicit none

    !global variables

    double precision,dimension(3)                           :: elps_tm
    double precision,dimension(3)                           :: t1s,t2s

    !local variables

    double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: u1r_rnd,u1i_rnd
    double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: u2r_rnd,u2i_rnd
    double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: u3r,u3i

    !timing variables
    double precision                :: elps_tm_mm
    double precision,dimension(2)   :: s_tm_mm, e_tm_mm
    double precision :: t1_mm,t2_mm

    double precision                :: elps_tm_u1, elps_tm_u2
    double precision,dimension(2)   :: s_tm_u1, e_tm_u1
    double precision,dimension(2)   :: s_tm_u2, e_tm_u2
    double precision :: t1_u1,t2_u1
    double precision :: t1_u2,t2_u2

    integer                         :: ic,jc,kc  !loop indices
    !start of the excution commands

    !initialising the fixed seed for the random number generator
    call init_fixed_random_seed(1234)

    !serial multiplication

    call cpu_time(t1_u1)
#if defined (LINUX)
    call gettimeofday_c(s_tm_u1)
#endif
    call su3random(u1r_rnd,u1i_rnd)
#if defined (LINUX)
    call gettimeofday_c(e_tm_u1)
    call elapsed_time_c(s_tm_u1,e_tm_u1,elps_tm_u1)
#endif
    call cpu_time(t2_u1)

    call cpu_time(t1_u2)
#if defined (LINUX)
    call gettimeofday_c(s_tm_u2)
#endif
    call su3random(u2r_rnd,u2i_rnd)
#if defined (LINUX)
    call gettimeofday_c(e_tm_u2)
    call elapsed_time_c(s_tm_u2,e_tm_u2,elps_tm_u2)
#endif
    call cpu_time(t2_u2)

    call cpu_time(t1_mm)
#if defined (LINUX)
    call gettimeofday_c(s_tm_mm)
#endif
    u3r = 0.0d0
    u3i = 0.0d0
    do ic=1,nc
       do jc=1,nc
          do kc=1,nc
             u3r(:,:,:,:,:,ic,jc) = u3r(:,:,:,:,:,ic,jc) + &
                  ( u1r_rnd(:,:,:,:,:,ic,kc) * u2r_rnd(:,:,:,:,:,kc,jc) - &
                  u1i_rnd(:,:,:,:,:,ic,kc) * u2i_rnd(:,:,:,:,:,kc,jc) )
             u3i(:,:,:,:,:,ic,jc) = u3i(:,:,:,:,:,ic,jc) + &
                  ( u1r_rnd(:,:,:,:,:,ic,kc) * u2i_rnd(:,:,:,:,:,kc,jc) + &
                  u1i_rnd(:,:,:,:,:,ic,kc) * u2r_rnd(:,:,:,:,:,kc,jc) )
          end do
       end do
    end do

#if defined (LINUX)
    call gettimeofday_c(e_tm_mm)
    call elapsed_time_c(s_tm_mm,e_tm_mm,elps_tm_mm)
#endif
    call cpu_time(t2_mm)

    elps_tm(1) = elps_tm_u1
    elps_tm(2) = elps_tm_u2
    elps_tm(3) = elps_tm_mm

    t1s(1) = t1_u1
    t1s(2) = t1_u2
    t1s(3) = t1_mm

    t2s(1) = t2_u1
    t2s(2) = t2_u2
    t2s(3) = t2_mm

    return
  end subroutine su3randomlinks_mul_timing

end module simple_SU3_tester
