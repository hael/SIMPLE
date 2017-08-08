! generic optimisation subroutines
module simple_opt_subs
implicit none

logical :: warn=.false.

contains

    !>  \brief  given an optimizer specification object spec, and n dimensional point spec%x
    !!          and an n dimensional direction spec%xi, linmin moves and resets
    !!          spec%x to where the function spec%costfun(spec%x) takes on a
    !!          minimum along the direction spec%xi from spec%x, and replaces
    !!          spec%xi by the actual vector displacement that spec%x was moved.
    !!          Returns as lowest_cost the value of spec%cosfun at the returned
    !!          location spec%x. All accomplished by calling the routines mnbrak
    !!          and brent. No derivatives needed
    subroutine linmin(spec,lowest_cost)
        use simple_opt_spec, only: opt_spec
        class(opt_spec), intent(inout) :: spec
        real,            intent(out)   :: lowest_cost
        integer :: j
        real    :: ax,bx,fa,fb,fx,xmin,xx
        if( spec%str_opt .ne. 'linmin' )then
            stop 'this optimzer specification object not 4 linmin; simple_opt_subs'
        endif
        if( .not. associated(spec%costfun) )then
            stop 'cost function not associated in opt_spec; linmin; simple_opt_subs'
        endif
        ax=0. ! initial guess for brackets
        xx=1.
        bx=2.
        call mnbrak(ax,xx,bx,fa,fx,fb)
        lowest_cost=brent(ax,xx,bx,xmin)
        do j=1,spec%ndim ! construct return vec
            spec%xi(j)=xmin*spec%xi(j)
            spec%x(j)=spec%x(j)+spec%xi(j)
        end do

        contains

            !>  \brief  routine for evaluating the move along the line direction
            function eval_move(x) result(r)
                real, intent(in) :: x
                real    :: r
                integer :: j
                do j=1,spec%ndim
                    spec%xt(j)=spec%x(j)+x*spec%xi(j)
                end do
                r = spec%costfun(spec%xt,spec%ndim)
                spec%nevals = spec%nevals+1
            end function eval_move

            !>  \brief  routine for initially bracketing a minimum
            !!          given the function eval_move(x), and given distinct
            !!          initial points ax and bx, this routine searches in the
            !!          downhill direction (defined by eval_move as evaluated at
            !!          the initial points) and returns new points ax, bx, cx
            !!          which bracket a minimum of the function. Function values
            !!          at the three points, fa, fb and fc are also returned
            !!          PARAMETERS: GOLD is the default ratio by which
            !!          successive intervals are magnified; GLIMIT is the
            !!          maximum magnification allowed for parabolic-fit step
            subroutine mnbrak(ax,bx,cx,fa,fb,fc)
                real, intent(inout) :: ax,bx,cx
                real, intent(out)   :: fa,fb,fc
                real, parameter     :: GOLD=1.618034, GLIMIT=100., TINY=1.e-10
                real :: dum,fu,q,r,u,ulim
                fa=eval_move(ax)
                fb=eval_move(bx)
                if(fb.gt.fa) then ! switch roles 4 a & b so that we can go downhill from a 2 b
                    dum=ax
                    ax=bx
                    bx=dum
                    dum=fb
                    fb=fa
                    fa=dum
                endif
                cx=bx+GOLD*(bx-ax) ! first guess 4 c
                fc=eval_move(cx)
                do while (fb.ge.fc)   ! keep returning here until we bracket
                    r=(bx-ax)*(fb-fc) ! compute u by parabolic extrapolation from a,b,c
                    q=(bx-cx)*(fb-fa) ! TINY prevents possible division by zero
                    u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
                    ulim=bx+GLIMIT*(cx-bx) ! we won't go farther than this. Test possibilities:
                    if((bx-u)*(u-cx).gt.0) then ! Parabolic u is between b and c: try it
                        fu=eval_move(u)
                        if(fu.lt.fc) then ! got a minimum between b and c
                            ax=bx
                            fa=fb
                            bx=u
                            fb=fu
                            cycle
                        else if(fu.gt.fb) then ! got a minimum between a and u
                            cx=u
                            fc=fu
                            cycle
                        endif
                        u=cx+GOLD*(cx-bx) ! parabolic fit was not useful. Use default magnification
                        fu=eval_move(u)
                    else if((cx-u)*(u-ulim).gt.0) then ! parabolic fit is between c and its allowed limit
                        fu=eval_move(u)
                        if(fu.lt.fc) then
                            bx=cx
                            cx=u
                            u=cx+GOLD*(cx-bx)
                            fb=fc
                            fc=fu
                            fu=eval_move(u)
                        endif
                    else if((u-ulim)*(ulim-cx).ge.0) then ! limit parabolic u to its max allowed val
                        u=ulim
                        fu=eval_move(u)
                    else ! reject parbolic u, use default magnification
                        u=cx+GOLD*(cx-bx)
                        fu=eval_move(u)
                    endif
                    ax=bx ! eliminate oldest point and continue
                    bx=cx
                    cx=u
                    fa=fb
                    fb=fc
                    fc=fu
                end do
            end subroutine mnbrak

            !>  \brief  Brent's method in one dimension
            !!          Given the eval_move function, and a bracketing triplet
            !!          of abscissas ax,bx,cx (such that bx is between ax and
            !!          cx, and f(bx) is less than both f(ax) and f(cx)), this
            !!          routine isolates the minimum to a fractional precision
            !!          of about spec%ftol using Brent's method. The abscissa of
            !!          the minimum is returned in xmin, and the minimum
            !!          function value is the returned function value
            function brent(ax,bx,cx,xmin) result(rbrent)
                use simple_math, only: shft
                real, intent(inout) :: ax,bx,cx
                real, intent(out)   :: xmin
                real, parameter :: CGOLD=.3819660, ZEPS=1e-10
                real :: rbrent,a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
                integer :: iter
                a=min(ax,cx) ! a & b must be in ascending order, though the input abscissas may not be
                b=max(ax,cx)
                v=bx         ! initializations
                w=v
                x=v
                e=0.         ! distance moved the step before last
                fx=eval_move(x)
                fv=fx
                fw=fx
                iter = 0
                do
                    iter = iter+1
                    if( iter == spec%maxits )then
!                        if(warn) write(*,'(a)') 'brent exceeding maximum iterations; simple_opt_subs'
                        return
                    endif
                    xm=0.5*(a+b)
                    tol1=spec%ftol*abs(x)+ZEPS
                    tol2=2.*tol1
                    if(abs(x-xm).le.(tol2-.5*(b-a)))then ! test for doneness here
                        xmin=x     ! exit with best values
                        rbrent=fx
                        return
                    endif
                    if(abs(e).gt.tol1)then ! construct a trial parabolic fit
                        r=(x-w)*(fx-fv)
                        q=(x-v)*(fx-fw)
                        p=(x-v)*q-(x-w)*r
                        q=.2*(q-r)
                        if(q.gt.0) p=-p
                        q=abs(q)
                        etemp=e
                        e=d
                        if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x))then ! acceptable parabolic fit?
                            e=merge(a-x,b-x,x >= xm)
                            d=CGOLD*e ! golden section step
                        else
                            d=p/q ! parabolic step
                            u=x+d
                            if(u-a.lt.tol2.or.b-u.lt.tol2) d=sign(tol1,xm-x)
                        endif
                    else
                        e=merge(a-x,b-x,x >= xm)
                        d=CGOLD*e
                    endif
                    u=merge(x+d,x+sign(tol1,d),abs(d) >= tol1)
                    fu = eval_move(u) ! this the one function evaluation per iteration
                    ! houskeeping follows...
                    if(fu.le.fx)then
                        if(u.ge.x)then
                            a=x
                        else
                            b=x
                        endif
                        call shft(v,w,x,u)
                        call shft(fv,fw,fx,fu)
                    else
                        if(u.lt.x)then
                            a=u
                        else
                            b=u
                        endif
                        if(fu.le.fw.or.w.eq.x)then
                            v=w
                            fv=fw
                            w=u
                            fw=fu
                        else if(fu.le.fv.or.v.eq.x.or.v.eq.w)then
                            v=u
                            fv=fu
                        endif
                    endif
                end do
            end function brent

    end subroutine linmin

    !>  \brief  line search routine (relying on gradient evals, used in BFGS)
    subroutine lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func,nevals)
        integer :: n, nevals
        logical :: check
        real    :: f, fold,stpmax,g(n),p(n)
        real    :: x(n),xold(n),ALF,TOLX
        parameter (ALF=1.e-4,TOLX=1.e-7)
        interface
            function func( vec, D ) result(cost)
                integer, intent(in) :: D
                real,    intent(in) :: vec(D)
                real                :: cost
            end function
        end interface
        integer :: i
        real    :: a,alam,alam2,alamin,b,disc,f2,fold2
        real    :: rhs1,rhs2,slope,sum,temp,test,tmplam
        check=.false.
        sum=0.;fold2=0.;alam=0.;alam2=0.;f2=0.
        do i=1,n
            sum=sum+p(i)*p(i)
        end do
        sum = sqrt(sum)
        if(sum.gt.stpmax)then ! scale if attempted step is too big
            do i=1,n
                p(i)=p(i)*stpmax/sum
            end do
        endif
        slope=0.
        do i=1,n
            slope=slope+g(i)*p(i)
        end do
        test=0. ! compute lambda_{min}
        do i=1,n
            temp=abs(p(i))/max(abs(xold(i)),1.)
            if(temp.gt.test)test=temp
        end do
        alamin=TOLX/test
        alam=1. ! always try full Newton step first
        do      ! start of loop
            do i=1,n
                x(i)=xold(i)+alam*p(i)
            end do
            f=func(x,n)
            nevals = nevals+1
            if(alam.lt.alamin)then ! convergence on delta_x. For zero finding,
                do i=1,n           ! the calling program should verify convergence
                    x(i)=xold(i)
                end do
                check=.true.
                return
            else if(f.le.fold+ALF*alam*slope)then ! sufficient function decrease
                return
            else ! backtrack
                if(alam.eq.1.)then ! first time
                    tmplam=-slope/(2.*(f-fold-slope))
                else               ! subsequent backtracks
                    rhs1=f-fold-alam*slope
                    rhs2=f2-fold2-alam2*slope
                    a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
                    b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
                    if(a.eq.0.)then
                        tmplam=-slope/(2.*b)
                    else
                        disc=b*b-3.*a*slope
                        if(disc.lt.0..and.warn) write(*,'(a)') 'roundoff problem in lnsrch; simple_opt_subs'
                        tmplam=(-b+sqrt(disc))/(3.*a)
                    endif
                    if(tmplam.gt..5*alam)tmplam=.5*alam
                endif
            endif
            alam2=alam
            f2=f
            fold2=fold
            alam=max(tmplam,.1*alam)
        end do
    end subroutine lnsrch

    !> \brief multidimensional minimization of the function func(x) (x(1:ndim)
    !>          is a vector in ndim dimensions) by the downhill simplex method
    !>          of Nelder and Mead.
    !!          The matrix p(1:ndim+1,1:ndim) is input/output. Its ndim+1 rows
    !!          are ndim-dimensional vectors which are the vertices of the
    !!          starting simplex. Input is also the vector y(1:ndim+1), whose
    !!          components must be pre-initialized to the values of funk
    !!          evaluated at the ndim+1 vertices (rows) of p. ftol is the
    !!          fractional convergence tolerance to be achieved in the function
    !!          value. On output, p and y will have been reset to ndim+1 new
    !!          points all within ftol of a minimum function value, and iter
    !!          gives the number of function evaluations taken.
    subroutine amoeba(p,y,pb,yb,ftol,func,iter,itmax,nevals)
        real,    intent(inout) :: p(:,:) !< the ndim+1 rows of p are ndim vec:s which are the vertices of the starting simplex
                                         !! the best point is put in slot 1 upon convergence
        real,    intent(inout) :: y(:)   !< must be pre-initialized to the values of func evaluated at the ndim+1 vertices (rows) of p
        real,    intent(inout) :: pb(:)  !< for updating the best point
        real,    intent(inout) :: yb     !< for updating the cost of best point
        real,    intent(in)    :: ftol   !< fractional convergence tolerance to be achieved in the function value (0.0005)
        integer, intent(out)   :: iter   !< number of exectuted iterations
        integer, intent(in)    :: itmax  !< maximum number of iterations
        integer, intent(inout) :: nevals !< number of costfun evals counter
        interface
            function func( vec, D ) result( cost ) !< the external function used for callback
                integer, intent(in) :: d
                real, intent(in)    :: vec(D)
                real :: cost
            end function
        end interface
        real, parameter :: tiny=1.0e-10
        integer         :: ihi,ndim
        real, dimension(size(p,2)) :: psum
        call amoeba_private

        contains

            subroutine amoeba_private
                 use simple_jiffys, only: assert_eq, swap
                integer :: i,ilo,inhi,loc(1)
                real :: rtol,ysave,ytry,ytmp
                ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,size(pb),'amoeba; simple_opt_subs')
                iter=0
                psum(:)=sum(p(:,:),dim=1)
                do
                    loc=minloc(y) ! determine which point has the highest (worst), next-highest, and lowest (best)
                    ilo=loc(1)
                    loc=maxloc(y)
                    ihi=loc(1)
                    ytmp=y(ihi)
                    y(ihi)=y(ilo)
                    loc=maxloc(y)
                    inhi=loc(1)
                    y(ihi)=ytmp
                    ! Compute the fractional range from highest to lowest and return if satisfactory
                    rtol=2.0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+tiny) ! relative tolerance
                    if(rtol < ftol)then
                        call swap(y(1),y(ilo))
                        call swap(p(1,:),p(ilo,:))
                        exit
                    end if
                    if(iter >= itmax)then
                        if(warn) write (*,*) 'itmax exceeded in amoeba; simple_opt_subs'
                        exit
                    endif
                    ! Begin a new iteration. First extrapolate by a factor −1 through the face of the
                    ! simplex across from the high point, i.e. reflect the simplex from the high point
                    ytry=amotry(-1.0)
                    iter=iter+1
                    if(ytry <= y(ilo))then
                        ! gives a result better than the best point, so try an additional extrapolation by a factor 2
                        ytry=amotry(2.0)
                        iter=iter+1
                    else if(ytry >= y(inhi))then
                        ! the reflected point is worse than the second-highest, so look for an intermediate
                        ! lower point, i.e. do a one-dimensional contraction
                        ysave=y(ihi)
                        ytry=amotry(0.5)
                        iter=iter+1
                        if(ytry >= ysave)then ! can’t seem to get rid of that high point
                            p(:,:)=0.5*(p(:,:)+spread(p(ilo,:),1,size(p,1))) ! better contract around the lowest (best) point
                            do i=1,ndim+1
                                if(i /= ilo)then
                                    y(i)=func(p(i,:),ndim)
                                    nevals = nevals+1
                                endif
                            end do
                            iter=iter+ndim ! keep track of function evaluations
                            psum(:)=sum(p(:,:),dim=1) ! recompute psum
                        end if
                    end if
                end do
                ! store best
                if( y(1) <= yb )then
                    pb(:) = p(1,:)
                    yb    = y(1)
                endif
            end subroutine amoeba_private

            !>  \brief  extrapolates by a factor fac through the face of the simplex across from the
            !!          high point, tries it, and replaces the high point if the new point is better
            function amotry(fac)
                real, intent(in) :: fac
                real :: amotry
                real :: fac1,fac2,ytry
                real, dimension(size(p,2)) :: ptry
                fac1=(1.0-fac)/ndim
                fac2=fac1-fac
                ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
                ytry=func(ptry,ndim)  ! evaluate the function at the trial point
                nevals = nevals+1
                if(ytry < y(ihi))then ! if it is better than the highest, then replace the highest
                    y(ihi)=ytry
                    psum(:)=psum(:)-p(ihi,:)+ptry(:)
                    p(ihi,:)=ptry(:)
                end if
                amotry=ytry
            end function amotry

    end subroutine amoeba

    !> \brief  stochastic hill climbing selection rule for distance arrays
    function shc_selector( dists, old_ind ) result( new_ind )
        use simple_ran_tabu, only: ran_tabu
        real, intent(in)    :: dists(:)
        integer, intent(in) :: old_ind
        type(ran_tabu)      :: rt
        integer :: i, n, new_ind, ntabu
        n     = size(dists)
        rt    = ran_tabu(n)
        ntabu = 0
        do i=1,n
            if( dists(i) < dists(old_ind) )then
                ! these are available (not tabu)
            else
                ! these are forbidden (tabu)
                call rt%insert(i)
                ntabu = ntabu+1
            endif
        end do
        if( ntabu == n )then
            new_ind = old_ind
        else
            new_ind = rt%irnd()
        endif
        call rt%kill
    end function shc_selector

    !> \brief  check the vector with respect to the limits
    subroutine check_and_correct_vec( spec, vec, corrected )
        use simple_opt_spec, only: opt_spec
        use simple_rnd,      only: ran3
        class(opt_spec), intent(in)    :: spec           !< specification
        real, intent(inout)            :: vec(spec%ndim) !< solution vector
        logical, intent(out), optional :: corrected   !< to indicate if vector was corrected
        integer                        :: j
        logical :: ccorrected
        ccorrected = .false.
        do j=1,spec%ndim
            if( spec%cyclic(j) ) then
                do while(vec(j) < spec%limits(j,1))
                    vec(j) = vec(j)+spec%limits(j,2)
                    ccorrected = .true.
                end do
                do while(vec(j) > spec%limits(j,2))
                    vec(j) = vec(j)-spec%limits(j,2)
                    ccorrected = .true.
                end do
            endif
        end do
        if( present(corrected) ) corrected = ccorrected
    end subroutine check_and_correct_vec

end module simple_opt_subs
