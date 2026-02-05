submodule (simple_oris) simple_oris_weights
use simple_ori_api
use simple_ori, only: ori
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine calc_hard_weights( self, frac )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: frac
        integer, allocatable :: order(:), states(:)
        integer :: i, j, lim, ind, n
        if( frac < 0.99 )then
            states = nint(self%get_all('state'))
            n      = count(states>0)
            lim    = nint(frac*real(n))
            order  = self%order()
            j = 0
            do i = 1,self%n
                ind = order(i)
                if( states(ind) == 0 )then
                    call self%o(ind)%set('w', 0.)
                else
                    j = j+1
                    if( j <= lim )then
                        call self%o(ind)%set('w', 1.)
                    else
                        call self%o(ind)%set('w', 0.)
                    endif
                endif
            end do
            deallocate(states,order)
        else
            call self%set_all2single('w', 1.)
        endif
    end subroutine calc_hard_weights

    module subroutine calc_soft_weights( self, frac )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: frac
        real,    allocatable :: scores(:), states(:), weights(:), weights_glob(:)
        integer, allocatable :: order(:), states_int(:)
        real    :: minw, mins
        integer :: i, lim, state, nstates
        if( self%isthere('corr') )then
            scores = self%get_all('corr')
            mins   = minval(scores, mask=scores > TINY)
            if( mins >= 0.8 )then
                call self%set_all2single('w', 1.0)
                return
            endif
            if( self%isthere('state') )then
                states = self%get_all('state')
            else
                allocate(states(self%n), source=1.0)
            endif
            nstates = nint(maxval(states))
            if( nstates == 1 )then
                weights = z_scores(scores, mask=scores > TINY .and. states > 0.5)
                minw    = minval(weights,  mask=scores > TINY .and. states > 0.5)
                where( scores > TINY .and. states > 0.5 )
                    weights = weights + abs(minw)
                elsewhere
                    weights = 0.
                endwhere
                call self%set_all('w', weights)
                deallocate(weights)
            else
                allocate(states_int(self%n),   source=nint(states))
                allocate(weights_glob(self%n), source=0.)
                do state=1,nstates
                    weights = z_scores(scores, mask=scores > TINY .and. states_int == state)
                    minw    = minval(weights,  mask=scores > TINY .and. states_int == state)
                    where( scores > TINY .and. states_int == state ) weights_glob = weights + abs(minw)
                    deallocate(weights)
                end do
                call self%set_all('w', weights_glob)
                deallocate(states_int, weights_glob)
            endif
            if( frac < 0.99 )then
                lim   = nint(frac*real(self%n))
                order = self%order()
                do i=1,self%n
                    if( i > lim ) call self%o(order(i))%set('w', 0.)
                end do
                deallocate(order)
            endif
            deallocate(scores, states)
        else
            call self%set_all2single('w', 1.0)
        endif
    end subroutine calc_soft_weights

    module subroutine calc_cavg_soft_weights( self, frac )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: frac
        real,    allocatable :: scores(:)
        integer, allocatable :: order(:), states(:)
        integer :: i, lim, state, nstates, nnonzero
        if( .not.self%isthere('corr') )then
            call self%set_all2single('w', 1.0)
            return
        endif
        scores = self%get_all('corr')
        if( self%isthere('state') )then
            states = self%get_all_asint('state')
            nstates = maxval(states)
        else
            allocate(states(self%n), source=1)
            nstates = 1
        endif
        if( minval(scores, mask=(scores>TINY).and.states>0) >= 0.85 )then
            call self%set_all2single('w', 1.0)
            return
        endif
        do state = 1,nstates
            call calc_weights( state )
        enddo
        if( frac < 0.999 )then
            order    = pack((/(i,i=1,self%n)/), mask=states > 0)
            nnonzero = size(order)
            lim      = nint((1.-frac)*real(nnonzero))
            scores   = scores(order(:))
            call hpsort(scores, order)
            do i = 1,lim,1
                call self%o(order(i))%set('w', 0.)
            end do
        endif
      contains

        subroutine calc_weights( s )
            integer,  intent(in) :: s
            real,    parameter   :: LOWER_BOUND_THRESHOLD = -4.0
            real,    allocatable :: weights(:), tmp(:)
            integer, allocatable :: inds(:), inds2(:)
            real    :: minw
            integer :: i, n
            n = count(states==s)
            if( n == 0 )then
                return
            else if( n < 5 )then
                inds = pack((/(i,i=1,self%n)/), mask=states==s)
                allocate(weights(n),source=1.0)
            else
                inds    = pack((/(i,i=1,self%n)/), mask=states==s)
                weights = scores(inds(:))
                inds2 = pack((/(i,i=1,n)/), mask=weights > TINY)
                tmp   = weights(inds2)
                where(weights <= TINY) weights = 0.
                call robust_scaling(tmp)
                minw = max(LOWER_BOUND_THRESHOLD,minval(tmp))
                tmp = (tmp - minw) / (1.0 - minw)
                do i = 1,size(inds2)
                    weights(inds2(i)) = min(1.0,max(0.0,tmp(i)))
                enddo
            endif
            do i = 1,n
                call self%o(inds(i))%set('w', weights(i))
            enddo
        end subroutine calc_weights

    end subroutine calc_cavg_soft_weights

    module subroutine calc_hard_weights2D( self, frac, ncls )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: frac
        integer,     intent(in)    :: ncls
        type(oris)           :: os
        type(ori)            :: o
        integer, allocatable :: pinds(:)
        integer :: i, icls, pop
        if( frac < 0.99 )then
            do icls=1,ncls
                call self%get_pinds(icls, 'class', pinds)
                if(.not.allocated(pinds)) cycle
                pop = size(pinds)
                call os%new(pop, self%o(1)%is_particle())
                do i=1,pop
                    call self%get_ori(pinds(i), o)
                    call os%set_ori(i, o)
                enddo
                call os%calc_hard_weights(frac)
                do i=1,pop
                    call self%set(pinds(i), 'w', os%get(i,'w'))
                enddo
                deallocate(pinds)
            enddo
        else
            call self%set_all2single('w', 1.)
        endif
        call o%kill
        call os%kill
    end subroutine calc_hard_weights2D

    module subroutine calc_soft_weights2D( self )
        class(oris), intent(inout) :: self
        real,    allocatable :: scores(:), states(:)
        real,    allocatable :: weights(:), weights_glob(:)
        integer, allocatable :: classes(:)
        integer :: icls, pop, ncls
        real    :: minw
        if( self%isthere('corr') )then
            scores = self%get_all('corr')
            classes    = nint(self%get_all('class'))
            if( self%isthere('state') )then
                states = self%get_all('state')
            else
                allocate(states(self%n), source=1.0)
            endif
            allocate(weights_glob(self%n), source=0.)
            ncls = maxval(classes)
            do icls=1,ncls
                pop = count(classes == icls .and. states > 0.5)
                if( pop == 0 )then
                    cycle
                else if( pop <= MINCLSPOPLIM )then
                    where( scores > TINY .and. (classes == icls .and. states > 0.5) ) weights_glob = 1.0
                else
                    if( count(scores > TINY .and. (classes == icls .and. states > 0.5)) <= MINCLSPOPLIM )then
                        where( scores > TINY .and. (classes == icls .and. states > 0.5) ) weights_glob = 1.0
                    else
                        weights = z_scores(scores, mask=scores > TINY .and. (classes == icls .and. states > 0.5))
                        minw    = minval(weights,      mask=scores > TINY .and. (classes == icls .and. states > 0.5))
                        where( scores > TINY .and. (classes == icls .and. states > 0.5) ) weights_glob = weights + abs(minw)
                        deallocate(weights)
                    endif
                endif
            end do
            call self%set_all('w', weights_glob)
            deallocate(scores, classes, states, weights_glob)
        else
            call self%set_all2single('w', 1.)
        endif
    end subroutine calc_soft_weights2D

    module subroutine find_best_classes( self, box, smpd, res_thresh,cls_mask, ndev )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: box
        real,        intent(in)    :: smpd, res_thresh, ndev
        logical,     intent(inout) :: cls_mask(1:self%n)
        real,    allocatable :: rfinds(:), corrs(:)
        logical, allocatable :: msk(:)
        real    :: ave, sdev, res, res_threshold, corr_threshold
        integer :: icls, nincl
        logical :: has_res, has_corr
        msk = cls_mask
        if( self%isthere('pop') )then
            do icls=1,self%n
                if(self%get(icls,'pop')<0.5) msk(icls) = .false.
            enddo
        endif
        has_corr = self%isthere('corr')
        if( has_corr )then
            do icls=1,self%n
                if(self%get(icls,'corr') < 0.0001) msk(icls) = .false.
            enddo
        endif
        has_res = self%isthere('res')
        if( has_res )then
            do icls=1,self%n
                if(self%get(icls,'res') > res_thresh) msk(icls) = .false.
            enddo
        endif
        nincl    = count(msk)
        cls_mask = msk
        if( has_res )then
            allocate(rfinds(self%n), source=0.)
            do icls=1,self%n
                res = self%get(icls, 'res')
                rfinds(icls) = real(calc_fourier_index(res,box,smpd))
            enddo
            ave  = sum(rfinds,mask=msk)/real(nincl)
            sdev = sqrt(sum((rfinds-ave)**2.,mask=msk)/real(nincl))
            res_threshold = max(ave-ndev*sdev,2.)
        else
            allocate(rfinds(self%n), source=huge(res_threshold))
            res_threshold = 0.
        endif
        if( has_corr )then
            corrs = self%get_all('corr')
            ave   = sum(corrs,mask=msk) / real(nincl)
            sdev  = sqrt(sum((corrs-ave)**2., mask=msk)/real(nincl))
            corr_threshold = ave-ndev*sdev
        else
            allocate(corrs(self%n), source=huge(corr_threshold))
            corr_threshold = 0.
        endif
        do icls=1,self%n
            if( cls_mask(icls) )then
                if(rfinds(icls)<res_threshold .and. corrs(icls)<corr_threshold)then
                    cls_mask(icls) = .false.
                endif
            endif
        enddo
        deallocate(msk,rfinds,corrs)
    end subroutine find_best_classes

end submodule simple_oris_weights
