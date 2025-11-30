program simple_test_complex_grad
! include 'simple_lib.f08'
! use simple_optimizer,         only: optimizer
! use simple_opt_factory,       only: opt_factory
! use simple_opt_spec,          only: opt_spec
! implicit none
! class(optimizer), pointer   :: opt_ptr=>null()      ! the generic optimizer object
! integer,          parameter :: NDIM=2, NRESTARTS=10, NPOINTS=10, NPI = 5
! type(opt_factory) :: ofac                           ! the optimization factory object
! type(opt_spec)    :: spec                           ! the optimizer specification object
! character(len=8)  :: str_opts                       ! string descriptors for the NOPTS optimizers
! real              :: lowest_cost, truth_xy(2), xy(2), lims(2,2), carg(NPOINTS), C(NPOINTS), T(NPOINTS),&
!                     &C_1st, T_1st, C_2nd, T_2nd, RHS_1st(NPI), RHS_2nd(NPI), x_comp, y_comp,&
!                     &cand1_x(NPI*NPI), cand1_y(NPI*NPI), cand2_x(NPI*NPI), cand2_y(NPI*NPI), minval, minx, miny, val
! complex           :: A(NPOINTS), B(NPOINTS), AB
! integer           :: i, j, cand1_cnt, cand2_cnt
! do i = 1, NPOINTS
!     A(i) = cmplx(ran3(), ran3()) * 5
!     C(i) = ran3() * PI
!     T(i) = ran3() * PI
! enddo
! truth_xy  = [ran3() * 4. - 2., ran3() * 4. - 2.]
! B         = A * cmplx(cos(C * truth_xy(1)), sin(T * truth_xy(2)))
! xy        = [ran3() - 0.5, ran3() - 0.5]
! str_opts  = 'lbfgsb'
! lims(1,1) = -5.
! lims(1,2) =  5.
! lims(2,1) = -5.
! lims(2,2) =  5.
! call spec%specify(str_opts, NDIM, limits=lims, nrestarts=NRESTARTS, factr  = 1.0d+5, pgtol = 1.0d-7)
! call spec%set_costfun(costfct)                                      ! set pointer to costfun
! call spec%set_gcostfun(gradfct)                                     ! set pointer to gradient of costfun
! call ofac%new(spec, opt_ptr)                                        ! generate optimizer object with the factory
! spec%x = xy
! call opt_ptr%minimize(spec, opt_ptr, lowest_cost)                   ! minimize the test function
! print *, 'CASE: quite stable complex-variable cost function'
! print *, 'starting xy = ', xy
! print *, 'truth    xy = ', truth_xy
! print *, 'searched xy = ', spec%x
! ! current shift searching cost function
! carg = C * truth_xy(1) + T * truth_xy(2)
! B    = A * cmplx(cos(carg), sin(carg))
! call spec%set_costfun(costfct_2)                                    ! set pointer to costfun
! call spec%set_gcostfun(gradfct_2)                                   ! set pointer to gradient of costfun
! spec%x = xy
! call opt_ptr%minimize(spec, opt_ptr, lowest_cost)                   ! minimize the test function
! print *, '-------------------------------------------------'
! print *, 'CASE: unstable complex-variable cost function'
! print *, 'starting xy = ', xy
! print *, 'truth    xy = ', truth_xy
! print *, 'searched xy = ', spec%x
! ! shift computation using two random equations
! RHS_1st = 0.
! C_1st   = 0.
! T_1st   = 0.
! i       = 1     ! first equation
! if( real(A(i)) < TINY .and. aimag(A(i)) < TINY )then
!     print *, 'bad choice of candidate for the equation'
!     stop
! endif
! AB      = B(i)/A(i)
! RHS_1st = RHS_1st + atan(aimag(AB)/real(AB))
! C_1st   = C_1st + C(i)
! T_1st   = T_1st + T(i)
! do j = 1, NPI
!     RHS_1st(j) = RHS_1st(j) + (j-1-NPI/2)*PI
! enddo
! RHS_2nd = 0.
! C_2nd   = 0.
! T_2nd   = 0.
! i       = 2     ! second equation
! if( real(A(i)) < TINY .and. aimag(A(i)) < TINY )then
!     print *, 'bad choice of candidate for the equation'
!     stop
! endif
! AB      = B(i)/A(i)
! RHS_2nd = RHS_2nd + atan(aimag(AB)/real(AB))
! C_2nd   = C_2nd + C(i)
! T_2nd   = T_2nd + T(i)
! do j = 1, NPI
!     RHS_2nd(j) = RHS_2nd(j) + (j-1-NPI/2)*PI
! enddo
! cand1_cnt = 0
! do i = 1, NPI
!     do j = 1, NPI
!         x_comp = -(RHS_1st(i)*T_2nd - RHS_2nd(j)*T_1st)/(T_1st*C_2nd - T_2nd*C_1st)
!         y_comp =  (RHS_1st(i)*C_2nd - RHS_2nd(j)*C_1st)/(T_1st*C_2nd - T_2nd*C_1st)
!         if( x_comp < -5. .or. x_comp > 5. .or. y_comp < -5. .or. y_comp > 5. )cycle
!         cand1_cnt = cand1_cnt + 1
!         cand1_x(cand1_cnt) = x_comp
!         cand1_y(cand1_cnt) = y_comp
!     enddo
! enddo
! ! another pair
! RHS_1st = 0.
! C_1st   = 0.
! T_1st   = 0.
! i       = 3     ! first equation
! if( real(A(i)) < TINY .and. aimag(A(i)) < TINY )then
!     print *, 'bad choice of candidate for the equation'
!     stop
! endif
! AB      = B(i)/A(i)
! RHS_1st = RHS_1st + atan(aimag(AB)/real(AB))
! C_1st   = C_1st + C(i)
! T_1st   = T_1st + T(i)
! do j = 1, NPI
!     RHS_1st(j) = RHS_1st(j) + (j-1-NPI/2)*PI
! enddo
! RHS_2nd = 0.
! C_2nd   = 0.
! T_2nd   = 0.
! i       = 4     ! second equation
! if( real(A(i)) < TINY .and. aimag(A(i)) < TINY )then
!     print *, 'bad choice of candidate for the equation'
!     stop
! endif
! AB      = B(i)/A(i)
! RHS_2nd = RHS_2nd(j) + atan(aimag(AB)/real(AB))
! C_2nd   = C_2nd + C(i)
! T_2nd   = T_2nd + T(i)
! do j = 1, NPI
!     RHS_2nd(j) = RHS_2nd(j) + (j-1-NPI/2)*PI
! enddo
! cand2_cnt = 0
! do i = 1, NPI
!     do j = 1, NPI
!         x_comp = -(RHS_1st(i)*T_2nd - RHS_2nd(j)*T_1st)/(T_1st*C_2nd - T_2nd*C_1st)
!         y_comp =  (RHS_1st(i)*C_2nd - RHS_2nd(j)*C_1st)/(T_1st*C_2nd - T_2nd*C_1st)
!         if( x_comp < -5. .or. x_comp > 5. .or. y_comp < -5. .or. y_comp > 5. )cycle
!         cand2_cnt = cand2_cnt + 1
!         cand2_x(cand2_cnt) = x_comp
!         cand2_y(cand2_cnt) = y_comp
!     enddo
! enddo
! ! finding the min between cand1 and cand2
! minval = huge(minval)
! do i = 1, cand1_cnt
!     do j = 1, cand2_cnt
!         val = sqrt((cand1_x(i)-cand2_x(j))**2 + (cand1_y(i)-cand2_y(j))**2)
!         if( val < minval )then
!             minval = val
!             minx   = (cand1_x(i)+cand2_x(j)) / 2.
!             miny   = (cand1_y(i)+cand2_y(j)) / 2.
!         endif
!     enddo
! enddo
! print *, '-------------------------------------------------'
! print *, 'minx = ', minx
! print *, 'miny = ', miny

! contains

!     function costfct( fun_self, x, d ) result( r )
!         class(*), intent(inout) :: fun_self
!         integer,  intent(in)    :: d
!         real,     intent(in)    :: x(d)
!         real    :: r
!         complex :: diff(NPOINTS)
!         diff = A * cmplx(cos(C * x(1)), sin(T * x(2))) - B
!         r    = sum(real(diff * conjg(diff)))
!     end function

!     subroutine gradfct( fun_self, x, grad, d )
!         class(*), intent(inout) :: fun_self
!         integer,  intent(in)    :: d
!         real,     intent(inout) :: x(d)
!         real,     intent(out)   :: grad(d)
!         complex :: diff(NPOINTS)
!         real    :: r
!         diff    = A * cmplx(cos(C * x(1)), sin(T * x(2))) - B
!         r       = sum(real(diff * conjg(diff)))
!         grad(1) = - 2. * sum(C * sin(C*x(1)) *  real(A * conjg(diff)))
!         grad(2) = - 2. * sum(T * cos(T*x(2)) * aimag(A * conjg(diff)))
!     end subroutine

!     function costfct_2( fun_self, x, d ) result( r )
!         class(*), intent(inout) :: fun_self
!         integer,  intent(in)    :: d
!         real,     intent(in)    :: x(d)
!         real    :: r, carg(NPOINTS)
!         complex :: diff(NPOINTS)
!         carg = C * x(1) + T * x(2)
!         diff = A * cmplx(cos(carg), sin(carg)) - B
!         r    = sum(real(diff * conjg(diff)))
!     end function

!     subroutine gradfct_2( fun_self, x, grad, d )
!         class(*), intent(inout) :: fun_self
!         integer,  intent(in)    :: d
!         real,     intent(inout) :: x(d)
!         real,     intent(out)   :: grad(d)
!         complex :: diff(NPOINTS)
!         real    :: r, carg(NPOINTS), tmp(NPOINTS)
!         carg    = C * x(1) + T * x(2)
!         diff    = A * cmplx(cos(carg), sin(carg)) - B
!         r       = sum(real(diff * conjg(diff)))
!         tmp     = -2. * (real(A * conjg(diff)) * sin(carg) + aimag(A * conjg(diff)) * cos(carg))
!         grad(1) = sum(C * tmp)
!         grad(2) = sum(T * tmp)
!     end subroutine

end program simple_test_complex_grad
