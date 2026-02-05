!@descr: routines for calculating statistics of oris (e.g. mean, median, min, max, etc.)
submodule (simple_oris) simple_oris_stats
use simple_ori_api
use simple_ori, only: ori
implicit none
#include "simple_local_flags.inc"

contains

    module function median_1( self, which ) result( med )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real,    allocatable :: vals(:), vals4med(:)
        logical, allocatable :: incl(:)
        real :: med
        incl     = self%included()
        vals     = self%get_all(which)
        vals4med = pack(vals, incl)
        med      = median_nocopy(vals4med)
    end function median_1

    module subroutine stats_1( self, which, ave, sdev, var, err )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: ave, sdev, var
        logical,          intent(out)   :: err
        real,    allocatable :: vals(:)
        integer, allocatable :: states(:)
        states        = self%get_all_asint('state')
        vals          = self%get_all(which)
        vals          = pack(vals, mask=states > 0)
        call moment(vals, ave, sdev, var, err)
        deallocate(vals, states)
    end subroutine stats_1

    module subroutine stats_2( self, which, statvars, mask, nozero )
        class(oris),        intent(inout) :: self
        character(len=*),   intent(in)    :: which
        type(stats_struct), intent(inout) :: statvars
        logical,            intent(in)    :: mask(self%n)
        logical, optional,  intent(in)    :: nozero
        real,    allocatable :: vals(:)
        integer, allocatable :: states(:)
        logical :: err, nnozero
        real    :: var
        nnozero = .false.
        if( present(nozero) ) nnozero = nozero
        vals          = self%get_all(which)
        states        = self%get_all_asint('state')
        vals          = pack(vals, states > 0 .and. mask)
        if( nnozero ) vals = pack(vals, mask=vals(:) > TINY)
        call moment(vals, statvars%avg, statvars%sdev, var, err)
        statvars%minv = minval(vals)
        statvars%maxv = maxval(vals)
        statvars%med  = median_nocopy(vals)
        deallocate(vals)
    end subroutine stats_2

    module subroutine stats_3( self, which, ave, sdev, var, err, fromto )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: ave, sdev, var
        logical,          intent(out)   :: err
        integer,          intent(in)    :: fromto(2)
        real,    allocatable            :: vals(:)
        integer, allocatable            :: states(:)
        states        = self%get_all_asint('state', fromto=fromto)
        vals          = self%get_all(which, fromto=fromto)
        vals          = pack(vals, mask=states > 0)
        call moment(vals, ave, sdev, var, err)
        deallocate(vals, states)
    end subroutine stats_3

    module subroutine minmax( self, which, minv, maxv )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: minv, maxv
        real    :: val, x
        integer :: i, mystate
        minv = huge(x)
        maxv = -huge(x)
        do i=1,self%n
            mystate = self%o(i)%get_state()
            if( mystate /= 0 )then
                val = self%o(i)%get(which)
                if( val < minv ) minv = val
                if( val > maxv ) maxv = val
            endif
        end do
    end subroutine minmax

    module subroutine spiral_1( self )
        class(oris), intent(inout) :: self
        real    :: h, theta, psi
        integer :: k
        if( self%n == 1 )then
            call self%o(1)%set_euler([0.,0.,0.])
        else if( self%n > 1 )then
            do k=1,self%n
                h     = -1.+2.*real(k-1)/real(self%n-1)
                theta = acos(h)
                if( k == 1 .or. k == self%n )then
                    psi = 0.
                else
                    psi = psi+3.6/(sqrt(real(self%n))*sqrt(1.-real(h)**2.))
                endif
                do while( psi > 2.*pi )
                    psi = psi-2.*pi
                end do
                call self%o(k)%set_euler([rad2deg(psi),rad2deg(theta),0.])
            end do
        else
            THROW_HARD('object nonexistent; spiral_1')
        endif
    end subroutine spiral_1

    module subroutine spiral_2( self, nsym, eullims )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nsym
        real,        intent(in)    :: eullims(3,2)
        logical, allocatable :: avail(:)
        type(oris) :: tmp
        integer    :: cnt, i, n, nprojs, lim
        real       :: e1lim, e2lim, frac, frac1, frac2
        if( nsym == 1 )then
            call self%spiral_1
            return
        endif
        e1lim  = eullims(1,2)
        e2lim  = eullims(2,2)
        frac1  = 360./e1lim
        frac2  = 1. / ( (1.-cos(deg2rad(e2lim))) /2. )
        frac   = frac1 * frac2
        n      = ceiling(real(self%n) * frac)
        call gen_c1
        nprojs = count(avail)
        if( nprojs < self%n )then
            n = n + self%n/2
            call gen_c1
            nprojs = count(avail)
        endif
        if( nprojs > self%n )then
            lim = max(2, floor(real(nprojs)/real(nprojs-self%n)))
            cnt  = 0
            do i = 1, n
                if(.not.avail(i))cycle
                cnt  = cnt + 1
                if(cnt == lim) then
                    avail(i) = .false.
                    cnt      = 0
                    nprojs   = nprojs-1
                    if( nprojs == self%n )exit
                endif
            enddo
        endif
        cnt = 0
        do i = 1, n
            if( avail(i) )then
                cnt = cnt + 1
                if(cnt > self%n)exit
                self%o(cnt) = tmp%o(i)
            endif
        enddo
        deallocate(avail)

        contains

            subroutine gen_c1
                integer :: i
                if( allocated(avail) )deallocate(avail)
                allocate(avail(n), source=.false.)
                call tmp%new(n, self%o(1)%is_particle())
                call tmp%spiral_1
                do i = 1, n
                    if( tmp%o(i)%e1get() <= e1lim .and. tmp%o(i)%e2get() <= e2lim )&
                    &avail(i) = .true.
                end do
            end subroutine gen_c1

    end subroutine spiral_2

    module function order( self ) result( inds )
        class(oris), intent(inout) :: self
        real,    allocatable :: scores(:)
        integer, allocatable :: inds(:)
        integer :: i
        inds   = (/(i,i=1,self%n)/)
        scores = self%get_all('corr')
        call hpsort(scores, inds)
        call reverse(inds)
        deallocate( scores )
    end function order

    module function order_cls( self, ncls ) result( inds )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: ncls
        integer, allocatable :: inds(:)
        real    :: classpops(ncls)
        integer :: i
        if(ncls <= 0) THROW_HARD('invalid number of classes; order_cls')
        allocate(inds(ncls))
        classpops = 0.0
        do i=1,ncls
            classpops(i) = real(self%get_pop(i, 'class'))
        end do
        inds = (/(i,i=1,ncls)/)
        call hpsort(classpops, inds)
        call reverse(inds)
    end function order_cls

    module subroutine detect_peaks( self, nnmat, corrs, peaks )
        class(oris), intent(in)    :: self
        integer,     intent(in)    :: nnmat(4,self%n)
        real,        intent(in)    :: corrs(self%n)
        logical,     intent(inout) :: peaks(self%n)
        real, allocatable :: corrs_packed(:)
        integer :: i, npeaks
        real    :: corr_t
        corr_t = 0.0
        do i = 1,self%n
            if( i /= nnmat(1,i) ) THROW_HARD('self is not set to the first entry of the 2nd dimension')
        end do
        do i = 1,self%n
            if( corrs(nnmat(1,i)) > corr_t )then
                peaks(i) = ( corrs(nnmat(1,i)) > corrs(nnmat(2,i)) .and.&
                &corrs(nnmat(1,i)) > corrs(nnmat(3,i)) .and.&
                &corrs(nnmat(1,i)) > corrs(nnmat(4,i)) )
            else
                peaks(i) = .false.
            endif
        end do
        npeaks = count(peaks .and. corrs > 0.)
        if( npeaks > 0 )then
            corrs_packed = pack(corrs, mask=peaks .and. corrs > 0.)
            call otsu(npeaks, corrs_packed, corr_t)
            deallocate(corrs_packed)
            where( corrs <= corr_t ) peaks = .false.
        endif
    end subroutine detect_peaks

    module function extremal_bound( self, thresh ) result( score_bound )
        class(oris),             intent(inout) :: self
        real,                    intent(in)    :: thresh
        real,    allocatable  :: scores(:), scores_incl(:)
        logical, allocatable  :: incl(:)
        integer :: n_incl, thresh_ind
        real    :: score_bound
        if( .not.self%isthere('corr') )then
            THROW_HARD('Metric: corr is unpopulated; extremal_bound')
        endif
        scores      = self%get_all('corr')
        incl        = self%included()
        scores_incl = pack(scores, mask=incl)
        n_incl      = size(scores_incl)
        call hpsort(scores_incl)
        thresh_ind  = nint(real(n_incl) * thresh)
        score_bound = scores_incl(thresh_ind)
        deallocate(scores, incl, scores_incl)
    end function extremal_bound

    module subroutine set_extremal_vars( self, extr_init, extr_iter, iter, frac_srch_space, do_extr, iextr_lim, update_frac )
        class(oris),       intent(in)  :: self
        real,              intent(in)  :: extr_init
        integer,           intent(in)  :: extr_iter, iter
        real,              intent(in)  :: frac_srch_space
        logical,           intent(out) :: do_extr
        integer,           intent(out) :: iextr_lim
        real,    optional, intent(in)  :: update_frac
        integer :: zero_pop
        logical :: l_update_frac
        do_extr           = .false.
        l_update_frac     = .false.
        if( present(update_frac) )then
            if( update_frac > 0.001 ) l_update_frac = .true.
        endif
        zero_pop  = self%n - self%get_noris(consider_state=.true.)
        iextr_lim = ceiling(2.*log(real(self%n-zero_pop)) * (2.-update_frac))
        if( l_update_frac )then
            iextr_lim = ceiling(2.*log(real(self%n-zero_pop)) * (2.-update_frac))
            if(iter==1 .or.(frac_srch_space <= 99. .and. extr_iter <= iextr_lim)) do_extr = .true.
        else
            iextr_lim = ceiling(2.*log(real(self%n-zero_pop)))
            if(iter==1 .or.(frac_srch_space <= 98. .and. extr_iter <= iextr_lim)) do_extr = .true.
        endif
    end subroutine set_extremal_vars

end submodule simple_oris_stats
