submodule (simple_oris) simple_oris_setters
use simple_ori_api
use simple_ori, only: ori
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine append_1( self, i, ori_in )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: ori_in
        integer,     intent(in)    :: i
        if( i < 0 .or. i > self%n )then
            THROW_WARN('index out of range; simple_oris % append')
            return
        endif
        call self%o(i)%append_ori(ori_in)
    end subroutine append_1

    module subroutine append_2( self1, self2 )
        class(oris), intent(inout) :: self1
        class(oris), intent(in)    :: self2
        integer :: i,nprev
        logical :: self1_isptcl, self2_isptcl
        if( self2%n == 0 ) return
        self2_isptcl = self2%o(1)%is_particle()
        if( self1%n == 0 )then
            call self1%copy(self2, self2_isptcl)
        else
            self1_isptcl = self1%o(1)%is_particle()
            if( self1_isptcl.eqv.self2_isptcl )then
                nprev = self1%n
                call self1%reallocate(self1%n+self2%n)
                do i = 1,self2%n
                    self1%o(nprev+i) = self2%o(i)
                end do
            else
                THROW_HARD('self1 and self2 do not have equivalent is_ptcl status')
            endif
        endif
    end subroutine append_2

    module subroutine copy_1( self_out, self_in, is_ptcl )
        class(oris), intent(inout) :: self_out
        class(oris), intent(in)    :: self_in
        logical,     intent(in)    :: is_ptcl
        integer :: i
        call self_out%new(self_in%n, is_ptcl)
        do i=1,self_in%n
            self_out%o(i) = self_in%o(i)
        end do
    end subroutine copy_1

    module subroutine copy_2( self_out, self_in )
        class(oris), intent(inout) :: self_out
        class(oris), intent(in)    :: self_in
        logical :: is_ptcl
        integer :: i
        if(self_in%get_noris() == 0) then
            call self_out%kill()
        else
            is_ptcl = self_in%is_particle()
            call self_out%new(self_in%n, is_ptcl)
            do i=1,self_in%n
                self_out%o(i) = self_in%o(i)
            end do
        end if
    end subroutine copy_2

    module subroutine reject( self, i )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        call self%o(i)%reject
    end subroutine reject

    module subroutine delete_entry_1( self, key )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer :: i
        do i=1,self%n
            call self%o(i)%delete_entry(key)
        end do
    end subroutine delete_entry_1

    module subroutine delete_entry_2( self, ind, key )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: key
        if( ind < 0 .or. ind > self%n )then
            THROW_WARN('index out of range; simple_oris % delete_entry_2')
            return
        endif
        call self%o(ind)%delete_entry(key)
    end subroutine delete_entry_2

    module subroutine delete_2Dclustering_1( self, keepshifts, keepcls )
        class(oris),       intent(inout) :: self
        logical, optional, intent(in)    :: keepshifts, keepcls
        integer :: i
        do i=1,self%n
            call self%o(i)%delete_2Dclustering(keepshifts, keepcls)
        end do
    end subroutine delete_2Dclustering_1

    module subroutine delete_2Dclustering_2( self, i, keepshifts, keepcls )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: i
        logical, optional, intent(in)    :: keepshifts, keepcls
        if( i < 0 .or. i > self%n )then
            THROW_WARN('index out of range; simple_oris % delete_2Dclustering_2')
            return
        endif
        call self%o(i)%delete_2Dclustering(keepshifts, keepcls)
    end subroutine delete_2Dclustering_2

    module subroutine delete_3Dalignment( self, keepshifts )
        class(oris),       intent(inout) :: self
        logical, optional, intent(in)    :: keepshifts
        integer :: i
        do i=1,self%n
            call self%o(i)%delete_3Dalignment(keepshifts)
        end do
    end subroutine delete_3Dalignment

    module subroutine delete( self, ind )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: ind
        type(ori),        allocatable   :: tmporis(:)
        if( ind < 0 .or. ind > self%n )then
            THROW_WARN('index out of range; simple_oris % delete')
            return
        endif
        allocate( tmporis(self%n - 1) )
        tmporis=[self%o(1:ind-1), self%o(ind+1:self%n)]
        call self%new(self%n - 1, .false.)
        self%o=tmporis
        if(allocated(tmporis)) deallocate(tmporis)
    end subroutine delete

    module subroutine transfer_2Dshifts( self_out, self_in )
        class(oris), intent(inout) :: self_out
        type(oris),   intent(in)   :: self_in
        integer :: i
        do i = 1,self_out%n
            call self_out%o(i)%set_shift(self_in%o(i)%get_2Dshift())
        enddo
    end subroutine transfer_2Dshifts

    module subroutine transfer_2Dparams_1( self, i, o_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        type(ori),   intent(in)    :: o_in
        call self%o(i)%transfer_2Dparams(o_in)
    end subroutine transfer_2Dparams_1

    module subroutine transfer_2Dparams_2( self, i, os_in, i_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, i_in
        class(oris), intent(in)    :: os_in
        call self%o(i)%transfer_2Dparams(os_in%o(i_in))
    end subroutine transfer_2Dparams_2

    module subroutine transfer_3Dparams_1( self, i, o_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        type(ori),   intent(in)    :: o_in
        call self%o(i)%transfer_3Dparams(o_in)
    end subroutine transfer_3Dparams_1

    module subroutine transfer_3Dparams_2( self, i, os_in, i_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, i_in
        class(oris), intent(in)    :: os_in
        call self%o(i)%transfer_3Dparams(os_in%o(i_in))
    end subroutine transfer_3Dparams_2

    module subroutine transfer_class_assignment( self_from, self_to )
        class(oris), intent(in)    :: self_from
        class(oris), intent(inout) :: self_to
        integer :: i
        if( self_from%n /= self_to%n ) THROW_HARD('Incongruent object instances')
        do i = 1,self_from%n
            call self_to%o(i)%set('class', self_from%o(i)%get_class())
        end do
    end subroutine transfer_class_assignment

    module subroutine set_euler( self, i, euls )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: euls(3)
        call self%o(i)%set_euler(euls)
    end subroutine set_euler

    module subroutine set_shift( self, i, vec )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: vec(2)
        call self%o(i)%set_shift(vec)
    end subroutine set_shift

    module subroutine set_state( self, i, state )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, state
        call self%o(i)%set_state(state)
    end subroutine set_state

    module subroutine set_class( self, i, cls )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, cls
        call self%o(i)%set_class(cls)
    end subroutine set_class

    module subroutine set_stkind( self, i, stkind )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, stkind
        call self%o(i)%set_stkind(stkind)
    end subroutine set_stkind

    module subroutine set_ogid( self, i, ogid )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, ogid
        call self%o(i)%set_ogid(ogid)
    end subroutine set_ogid

    module subroutine e1set( self, i, e1 )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: e1
        call self%o(i)%e1set(e1)
    end subroutine e1set

    module subroutine e2set( self, i, e2 )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: e2
        call self%o(i)%e2set(e2)
    end subroutine e2set

    module subroutine e3set( self, i, e3 )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: e3
        call self%o(i)%e3set(e3)
    end subroutine e3set

    module subroutine set_1( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real,             intent(in)    :: val
        call self%o(i)%set(key, val)
    end subroutine set_1

    module subroutine set_2( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        character(len=*), intent(in)    :: val
        call self%o(i)%set(key, val)
    end subroutine set_2

    module subroutine set_3( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        integer,          intent(in)    :: val
        call self%o(i)%set(key, val)
    end subroutine set_3

    module subroutine set_4( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real(dp),         intent(in)    :: val
        call self%o(i)%set(key, val)
    end subroutine set_4
    
    module subroutine set_5( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        class(string),    intent(in)    :: val
        call self%o(i)%set(key, val)
    end subroutine set_5

    module subroutine set_dfx( self, i, dfx )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: dfx
        call self%o(i)%set_dfx(dfx)
    end subroutine set_dfx

    module subroutine set_dfy( self, i, dfy )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: dfy
        call self%o(i)%set_dfy(dfy)
    end subroutine set_dfy

    module subroutine set_ori( self, i, o )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        class(ori),  intent(in)    :: o
        self%o(i) = o
    end subroutine set_ori

    module subroutine transfer_ori( self, i, self2transfer, i2transfer )
        class(oris), intent(inout) :: self
        class(oris), intent(in)    :: self2transfer
        integer,     intent(in)    :: i, i2transfer
        self%o(i) = self2transfer%o(i2transfer)
    end subroutine transfer_ori

    module subroutine set_all_1( self, which, vals )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: vals(self%n)
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, vals(i))
        enddo
    end subroutine set_all_1

    module subroutine set_all_2( self, which, vals )
        class(oris),           intent(inout) :: self
        character(len=*),      intent(in)    :: which
        character(len=STDLEN), intent(in)    :: vals(self%n)
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, vals(i))
        enddo
    end subroutine set_all_2

    module subroutine set_all_3( self, which, vals )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: vals(self%n)
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, vals(i))
        enddo
    end subroutine set_all_3

    module subroutine set_all2single_1( self, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: val
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, val)
        end do
    end subroutine set_all2single_1

    module subroutine set_all2single_2( self, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        character(len=*), intent(in)    :: val
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, val)
        end do
    end subroutine set_all2single_2

    module subroutine set_all2single_3( self, which, ival )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: ival
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, ival)
        end do
    end subroutine set_all2single_3

    module subroutine set_field2single_1( self, field, ind, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: field
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: val
        integer :: i
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i=1,self%n
            if( self%o(i)%get(trim(field)) == ind ) call self%o(i)%set(which, val)
        end do
        !$omp end parallel do 
    end subroutine set_field2single_1

    module subroutine set_field2single_2( self, field, ind, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: field
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: which
        character(len=*), intent(in)    :: val
        integer :: i
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i=1,self%n
            if( self%o(i)%get(trim(field)) == ind ) call self%o(i)%set(which, val)
        end do
        !$omp end parallel do 
    end subroutine set_field2single_2

    module subroutine set_field2single_3( self, field, ind, which, ival )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: field
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: ival
        integer :: i
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i=1,self%n
            if( self%o(i)%get(trim(field)) == ind ) call self%o(i)%set(which, ival)
        end do
        !$omp end parallel do 
    end subroutine set_field2single_3

    module subroutine set_projs( self, e_space )
        class(oris), intent(inout) :: self
        class(oris), intent(in)    :: e_space
        integer :: i
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i=1,self%n
            call self%set(i, 'proj', e_space%find_closest_proj(self%o(i)))
        end do
        !$omp end parallel do
    end subroutine set_projs

    module subroutine remap_projs( self, e_space, mapped_projs )
        class(oris), intent(in)  :: self
        class(oris), intent(in)  :: e_space
        integer,     intent(out) :: mapped_projs(self%n)
        integer :: i
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i=1,self%n
            mapped_projs(i) = e_space%find_closest_proj(self%o(i))
        end do
        !$omp end parallel do
    end subroutine remap_projs

    module subroutine proj2class( self )
        class(oris), intent(inout) :: self
        integer :: i
        if( .not. self%isthere('proj') ) THROW_HARD('No proj indices to turn into class indices; proj2class')
        do i=1,self%n
            call self%o(i)%set('class', self%o(i)%get('proj'))
        end do
    end subroutine proj2class

    module subroutine e3swapsgn( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%e3set(i,360.-self%e3get(i))
        end do
    end subroutine e3swapsgn

    module subroutine swape1e3( self )
        class(oris), intent(inout) :: self
        integer :: i
        real :: e, euls(3)
        do i=1,self%n
            euls = self%get_euler(i)
            e = euls(1)
            euls(1) = euls(3)
            euls(3) = e
            call self%set_euler(i,euls)
        end do
    end subroutine swape1e3

    module subroutine zero( self, which )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, 0.)
        end do
    end subroutine zero

    module subroutine zero_projs( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%e1set(0.)
            call self%o(i)%e2set(0.)
        end do
    end subroutine zero_projs

    module subroutine zero_shifts( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%set('x', 0.)
            call self%o(i)%set('y', 0.)
        end do
    end subroutine zero_shifts

    module subroutine zero_inpl( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%e3set(0.)
            call self%o(i)%set('x', 0.)
            call self%o(i)%set('y', 0.)
        end do
    end subroutine zero_inpl

    module subroutine mul_shifts( self, mul )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: mul
        integer :: i
        do i=1,self%n
            call self%o(i)%set_shift(mul*self%o(i)%get_2Dshift())
        end do
    end subroutine mul_shifts

    module subroutine rnd_oris( self, trs, eullims )
        class(oris),    intent(inout) :: self
        real, optional, intent(in)    :: trs
        real, optional, intent(inout) :: eullims(3,2)
        integer :: i
        do i=1,self%n
            call self%o(i)%rnd_ori(trs, eullims)
        end do
    end subroutine rnd_oris

    module subroutine gau_rnd_shifts( self, std )
        class(oris),    intent(inout) :: self
        real,           intent(in)    :: std
        integer :: i
        do i=1,self%n
            call self%o(i)%gau_rnd_shift(std)
        end do
    end subroutine gau_rnd_shifts

    module subroutine rnd_ori( self, i, trs, eullims )
        class(oris),    intent(inout) :: self
        integer,        intent(in)    :: i
        real, optional, intent(in)    :: trs
        real, optional, intent(inout) :: eullims(3,2)
        call self%o(i)%rnd_ori( trs, eullims )
    end subroutine rnd_ori

    module subroutine rnd_inpls( self, trs )
        class(oris),    intent(inout) :: self
        real, optional, intent(in)    :: trs
        integer :: i
        real :: x, y
        do i=1,self%n
            if( present(trs) )then
                if( abs(trs) < TINY )then
                    x = 0.
                    y = 0.
                else
                    x = ran3()*2.0*trs-trs
                    y = ran3()*2.0*trs-trs
                endif
            else
                x = 0.
                y = 0.
            endif
            call self%o(i)%e3set(ran3()*359.99)
            call self%o(i)%set('x', x)
            call self%o(i)%set('y', y)
        end do
    end subroutine rnd_inpls

    module subroutine rnd_ctf( self, kv, cs, fraca, defocus, deferr, astigerr )
        class(oris),    intent(inout) :: self
        real,           intent(in)    :: kv, cs, fraca, defocus, deferr
        real, optional, intent(in)    :: astigerr
        integer :: i
        real    :: dfx, dfy, angast, err
        do i=1,self%n
            call self%o(i)%set('kv',    kv   )
            call self%o(i)%set('cs',    cs   )
            call self%o(i)%set('fraca', fraca)
            do
                err = ran3()*deferr
                if( ran3() < 0.5 )then
                    dfx = defocus-err
                else
                    dfx = defocus+err
                endif
                if( dfx > 0. ) exit
            end do
            call self%o(i)%set_dfx(dfx)
            if( present(astigerr) )then
                do
                    err = ran3()*astigerr
                    if( ran3() < 0.5 )then
                        dfy = dfx-err
                    else
                        dfy = dfx+err
                    endif
                    if( dfy > 0. ) exit
                end do
                angast = ran3()*359.99
                call self%o(i)%set_dfy(dfy)
                call self%o(i)%set('angast', angast)
            endif
        end do
    end subroutine rnd_ctf

    module subroutine rnd_states( self, nstates )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nstates
        integer, allocatable       :: states(:)
        type(ran_tabu)             :: rt
        integer :: i, state
        if( nstates > 1 )then
            allocate( states(self%n) )
            rt = ran_tabu(self%n)
            call rt%balanced(nstates, states)
            do i=1,self%n
                state = self%o(i)%get_state()
                if( state /= 0 )then
                    call self%o(i)%set_state(states(i))
                endif
            end do
            call rt%kill
            deallocate(states)
        else if( nstates<=0)then
            THROW_HARD('invalid value for nstates; rnd_states')
        else
            ! nstates = 1; zero-preserving
            do i=1,self%n
                state = self%o(i)%get_state()
                if(  state /= 0 ) call self%o(i)%set_state(1)
            end do
        endif
    end subroutine rnd_states

    module subroutine rnd_lps( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%set('lp', ran3()*100.)
        end do
    end subroutine rnd_lps

    module subroutine rnd_corrs( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%set('corr', ran3())
        end do
    end subroutine rnd_corrs

    module subroutine rnd_proj_space( self, nsample, o_prev, thres, eullims )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: nsample
        class(ori), optional, intent(inout) :: o_prev
        real,       optional, intent(inout) :: eullims(3,2)
        real,       optional, intent(in)    :: thres
        type(ori) :: o_stoch
        integer   :: i
        logical   :: within_lims, found
        if( present(o_prev) .and. .not.present(thres) ) &
            & THROW_HARD('missing angular threshold in rnd_proj_space')
        if( .not.present(o_prev) .and. present(thres) ) &
            & THROW_HARD('missing orientation in rnd_proj_space')
        within_lims = .false.
        if( present(eullims) )within_lims = .true.
        call self%new(nsample, is_ptcl=.false.)
        if( present(o_prev).and.present(thres) )then
            do i=1,self%n
                found = .false.
                do while( .not.found )
                    call o_stoch%new(is_ptcl=.false.)
                    if( within_lims )then
                        call o_stoch%rnd_euler( eullims )
                    else
                        call o_stoch%rnd_euler
                    endif
                    if( rad2deg( o_stoch.euldist.o_prev ) > thres )cycle
                    found = .true.
                    call o_stoch%e3set( 0.)
                    self%o( i ) = o_stoch
                end do
            end do
        else
            do i=1,self%n
                if( within_lims)then
                    call self%rnd_ori( i, eullims=eullims )
                else
                    call self%rnd_ori( i )
                endif
                call self%e3set(i,0.)
            end do
        endif
    end subroutine rnd_proj_space

    module subroutine revshsgn( self )
        class(oris), intent(inout) :: self
        integer :: i
        real :: x, y
        do i=1,self%n
            x = self%o(i)%get('x')
            y = self%o(i)%get('y')
            call self%o(i)%set('x', -x)
            call self%o(i)%set('y', -y)
        end do
    end subroutine revshsgn

    module subroutine revorisgn( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%set('x', -self%o(i)%get('x'))
            call self%o(i)%set('y', -self%o(i)%get('y'))
            call self%o(i)%set_euler(-self%o(i)%get_euler())
        end do
    end subroutine revorisgn

    module subroutine ini_tseries( self, nsplit, state_or_class )
        use simple_map_reduce, only: split_nobjs_even
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: nsplit
        character(len=*), intent(in)    :: state_or_class
        integer, allocatable :: parts(:,:)
        integer :: ipart, iptcl
        parts = split_nobjs_even(self%n, nsplit)
        do ipart=1,nsplit
            do iptcl=parts(ipart,1),parts(ipart,2)
                call self%o(iptcl)%set(trim(state_or_class), ipart)
            end do
        end do
        if(allocated(parts))deallocate(parts)
    end subroutine ini_tseries

    module subroutine symmetrize( self, nsym )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nsym
        type(oris) :: tmp
        type(ori)  :: o
        integer    :: cnt, i, j
        logical    :: is_ptcl
        is_ptcl = self%o(1)%is_particle()
        tmp = oris(self%get_noris()*nsym, is_ptcl)
        cnt = 0
        do i=1,self%get_noris()
            do j=1,nsym
                cnt = cnt+1
                call self%get_ori(i, o)
                call tmp%set_ori(cnt, o)
            end do
        end do
        call self%copy(tmp, is_ptcl)
        call tmp%kill
        call o%kill
    end subroutine symmetrize

    module subroutine merge( self1, self2add )
        class(oris), intent(inout) :: self1, self2add
        type(oris) :: self
        integer    :: ntot, cnt, i
        logical    :: is_ptcl
        is_ptcl = self1%o(1)%is_particle()
        ntot    = self1%n+self2add%n
        self    = oris(ntot, is_ptcl)
        cnt     = 0
        if( self1%n > 0 )then
            do i=1,self1%n
                cnt = cnt+1
                self%o(cnt) = self1%o(i)
            end do
        endif
        if( self2add%n > 0 )then
            do i=1,self2add%n
                cnt = cnt+1
                self%o(cnt) = self2add%o(i)
            end do
        endif
        call self1%copy(self, is_ptcl)
        call self%kill
        call self2add%kill
    end subroutine merge

    module subroutine partition_eo( self )
        class(oris),       intent(inout) :: self
        integer :: i
        do i = 1, self%n-1, 2
            call self%set(i,  'eo', 0.)
            call self%set(i+1,'eo', 1.)
        end do
        if(is_even(self%n))then
            call self%set(self%n,'eo',1.)
        else
            call self%set(self%n,'eo',0.)
        endif
    end subroutine partition_eo

    module subroutine str2ori( self, i, line )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: line
        call self%o(i)%str2ori(line, self%o(1)%is_particle())
    end subroutine str2ori

    module subroutine str2ori_ctfparams_state_eo( self, i, line )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: line
        type(ori) :: o_tmp
        call o_tmp%str2ori(line, self%o(1)%is_particle())
        if( o_tmp%isthere('smpd')    ) call self%o(i)%set('smpd',    o_tmp%get('smpd'))
        if( o_tmp%isthere('kv')      ) call self%o(i)%set('kv',      o_tmp%get('kv'))
        if( o_tmp%isthere('cs')      ) call self%o(i)%set('cs',      o_tmp%get('cs'))
        if( o_tmp%isthere('fraca')   ) call self%o(i)%set('fraca',   o_tmp%get('fraca'))
        if( o_tmp%isthere('phshift') ) call self%o(i)%set('phshift', o_tmp%get('phshift'))
        if( o_tmp%isthere('dfx')     ) call self%o(i)%set_dfx(       o_tmp%get_dfx())
        if( o_tmp%isthere('dfy')     ) call self%o(i)%set_dfy(       o_tmp%get_dfy())
        if( o_tmp%isthere('angast')  ) call self%o(i)%set('angast',  o_tmp%get('angast'))
        if( o_tmp%isthere('state')   )then
            call self%o(i)%set('state', o_tmp%get('state'))
        else
            call self%o(i)%set('state', 1.0)
        endif
        if( o_tmp%isthere('eo') ) call self%o(i)%set('eo',      o_tmp%get('eo'))
        call o_tmp%kill
    end subroutine str2ori_ctfparams_state_eo

    module subroutine set_ctfvars( self, i, ctfvars )
        class(oris),     intent(inout) :: self
        integer,         intent(in)    :: i
        type(ctfparams), intent(in)    :: ctfvars
        call self%o(i)%set_ctfvars(ctfvars)
    end subroutine set_ctfvars

    module subroutine gen_balanced_partitions( self, nparts, parts, err )
        class(oris),          intent(in)    :: self
        integer,              intent(in)    :: nparts
        integer, allocatable, intent(inout) :: parts(:,:)
        logical,              intent(out)   :: err
        real, allocatable :: rstates(:)
        integer :: nobjs_per_part, i, ipart, m
        err = .false.
        if( allocated(parts) ) deallocate(parts)
        if( .not.self%isthere('state') )then
            err = .true.
            return
        endif
        allocate(parts(nparts,2),source=0)
        nobjs_per_part = ceiling(real(self%n)/real(nparts))
        rstates    = self%get_all('state')
        parts(1,1) = 1
        ipart      = 1
        m          = 0
        do i = 1,self%n
            if( rstates(i) < 0.5 )cycle
            m = m+1
            if( m == nobjs_per_part )then
                m = 0
                parts(ipart,2) = i
                if( ipart < nparts ) parts(ipart+1,1) = i+1
                ipart = ipart + 1
                if( ipart == nparts )exit
            endif
        enddo
        parts(nparts,2) = self%n
        deallocate(rstates)
    end subroutine gen_balanced_partitions

end submodule simple_oris_setters