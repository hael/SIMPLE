!@descr: the in-plane marginal of a probability table's batch for the students' teacher: each entry drawn as the fill draws it, its in-plane row weighed as the in-plane posterior insertion weighs it, and the in-plane posteriors of a row's heaviest references
! A table entry is the loss at one in-plane rotation drawn from its row's posterior at the table's temperature
! (pftc%gen_posterior_val); here the row is computed once and drawn with the same keyed variate (the same entry, bit
! for bit). The same row is weighed by the insertion's own quadrature (inpl_quadrature) at the insertion's physical
! temperature: its log partition sum, minimum and in-plane posterior are the insertion's for the shift it would take
! the entry at, and the shift prior's term of the view weight is the caller's (0 in 2D). The marginal mass
! of a reference is the insertion's view weight exp(lnz - T dmin - pen) (the class student's teacher,
! doc/algorithms/class_student.md); only the class student's learn program and its validation reference ask for it.
! audit_entry recomputes one entry by the insertion's route and reports the differences.
module simple_eul_prob_marginal
use simple_defs,               only: sp, dp, INPL_POSTERIOR_BUDGET
use simple_error,              only: simple_exception
use simple_polarft_calc,       only: polarft_calc, SHERRSQ
use simple_eul_prob_tab_utils, only: draw_posterior_row, inpl_quadrature, inpl_view_weights, prob_row_marginal, &
    &marginal_dist, posterior_temperature, temperature_or_argmin
use simple_rnd,                only: rnd_key, rnd_keyed, rnd_draw_key, RND_FAM_INPL
use simple_srch_sort_loc,      only: maxnloc
implicit none

public :: eul_prob_marginal, inpl_insert_temp
private
#include "simple_local_flags.inc"

type :: eul_prob_marginal
    private
    integer                  :: nrefs  = 0, nrots = 0, kin = 0, i_first = 0, i_last = -1
    real(dp),    allocatable :: lnz(:,:)             !< (nrefs, i_first:i_last) ln of the in-plane partition sum
    real(dp),    allocatable :: pen(:,:)             !< (nrefs, i_first:i_last) the view weight's shift-prior term
    real,        allocatable :: dmin(:,:)            !< (nrefs, i_first:i_last) its row's minimum
    real,        allocatable :: sh(:,:,:)            !< (2, nrefs, i_first:i_last) the shift its row was last computed at
    real,        allocatable :: temp(:)              !< (i_first:i_last) the physical temperature (0 for the argmin)
    integer,     allocatable :: nin(:)               !< (i_first:i_last) the in-plane rows kept
    integer,     allocatable :: refs(:,:)            !< (kin, i_first:i_last) their references, in descending marginal mass
    real,        allocatable :: rows(:,:,:)          !< (nrots, kin, i_first:i_last) their in-plane posteriors
    real,        allocatable :: wrow(:,:), wcdf(:,:) !< (nrots, nthr) per thread: a row, the draw's work
    complex(sp), allocatable :: bank(:,:,:)          !< (pftsz+1, nk, nthr) per thread: the shifted particle
    logical                  :: exists = .false.
contains
    procedure :: new
    procedure :: score_entry
    procedure :: rescore_entry
    procedure :: finish_row
    procedure :: get_marginal
    procedure :: audit_entry
    procedure :: audit_view_weights
    procedure :: kill
end type eul_prob_marginal

contains

    !> the in-plane posterior insertion's temperature of particle iptcl (insert_inplane): the physical
    !! one, the argmin's marker without a scale
    real function inpl_insert_temp( pftc, iptcl ) result( temp )
        class(polarft_calc), intent(in) :: pftc
        integer,             intent(in) :: iptcl
        temp = temperature_or_argmin(posterior_temperature(pftc%wsqsum_ptcl_now(iptcl), pftc%get_pftsz()))
    end function inpl_insert_temp

    !> the batch i_first..i_last of a table of nrefs references on pftc's rotation grid, the kin heaviest references'
    !! in-plane rows kept, nthr threads' work (the team the fill runs)
    subroutine new( self, pftc, nrefs, kin, i_first, i_last, nthr )
        class(eul_prob_marginal), intent(inout) :: self
        class(polarft_calc),      intent(in)    :: pftc
        integer,                  intent(in)    :: nrefs, kin, i_first, i_last, nthr
        integer :: kfromto(2)
        call self%kill
        if( nrefs < 1 .or. kin < 1 .or. i_last < i_first .or. nthr < 1 ) THROW_HARD('invalid batch; eul_prob_marginal%new')
        self%nrefs   = nrefs
        self%nrots   = pftc%get_nrots()
        self%kin     = min(kin, nrefs)
        self%i_first = i_first
        self%i_last  = i_last
        kfromto      = pftc%get_kfromto()
        allocate(self%lnz(nrefs, i_first:i_last), self%pen(nrefs, i_first:i_last), source=0.d0)
        allocate(self%dmin(nrefs, i_first:i_last), source=huge(1.))
        allocate(self%sh(2, nrefs, i_first:i_last), source=0.)
        allocate(self%temp(i_first:i_last), source=0.)
        allocate(self%nin(i_first:i_last), source=0)
        allocate(self%refs(self%kin, i_first:i_last), source=0)
        allocate(self%rows(self%nrots, self%kin, i_first:i_last), source=0.)
        allocate(self%wrow(self%nrots, nthr), self%wcdf(self%nrots, nthr), source=0.)
        allocate(self%bank(2*(kfromto(2)-kfromto(1)+1), pftc%get_pftsz()+1, nthr))
        self%exists = .true.
    end subroutine new

    !> Entry ri of row i: reference iref against particle iptcl shifted by shift (the shifted-reference kernel) or, with
    !! bank, the particle the caller shifted by shift once; the in-plane rotation drawn from the row's posterior at
    !! tdraw with the variate pftc%gen_posterior_val keys (scope, particle, reference), so dist and irot are its, bit
    !! for bit; the row's in-plane marginal at the insertion's temperature kept. On thread ithr's work
    subroutine score_entry( self, pftc, scope, iref, ri, i, iptcl, shift, tdraw, ithr, dist, irot, bank )
        class(eul_prob_marginal), intent(inout)           :: self
        class(polarft_calc),      intent(inout)           :: pftc
        type(rnd_key),            intent(in)              :: scope
        integer,                  intent(in)              :: iref, ri, i, iptcl, ithr
        real,                     intent(in)              :: shift(2), tdraw
        real,                     intent(out)             :: dist
        integer,                  intent(out)             :: irot
        complex(sp),              optional,    intent(in) :: bank(:,:)
        real(dp) :: u
        if( present(bank) )then
            call pftc%gen_euclid_dists_bank(iref, iptcl, bank, self%wrow(:,ithr))
        else
            call pftc%gen_euclid_dists(iref, iptcl, shift, self%wrow(:,ithr))
        endif
        u = rnd_keyed(rnd_draw_key(scope, iptcl, RND_FAM_INPL, reference=iref))
        call draw_posterior_row(self%nrots, self%wrow(:,ithr), tdraw, self%wcdf(:,ithr), u, irot, dist)
        if( irot < 1 )then
            irot = 1
            dist = huge(dist)
        endif
        self%sh(:,ri,i) = shift
        call weigh_row(self, ri, i, inpl_insert_temp(pftc, iptcl), ithr)
    end subroutine score_entry

    !> entry ri of row i after the fill refined it at shift: its in-plane row there (as the in-plane insertions compute
    !! it) replaces the marginal
    subroutine rescore_entry( self, pftc, iref, ri, i, iptcl, shift, ithr )
        class(eul_prob_marginal), intent(inout) :: self
        class(polarft_calc),      intent(inout) :: pftc
        integer,                  intent(in)    :: iref, ri, i, iptcl, ithr
        real,                     intent(in)    :: shift(2)
        call inpl_row_at(self, pftc, iref, iptcl, shift, ithr)
        self%sh(:,ri,i) = shift
        call weigh_row(self, ri, i, inpl_insert_temp(pftc, iptcl), ithr)
    end subroutine rescore_entry

    !> row i complete, pens(nrefs) the view weights' shift-prior terms: the kin references of largest marginal mass (the
    !! least marginal loss; a reference without a valid row excluded) and their in-plane posteriors at the insertion's
    !! temperature, each at the shift its entry was last scored at (iref_of(ri) the calculator's reference)
    subroutine finish_row( self, pftc, i, iptcl, iref_of, pens, ithr )
        class(eul_prob_marginal), intent(inout) :: self
        class(polarft_calc),      intent(inout) :: pftc
        integer,                  intent(in)    :: i, iptcl, iref_of(:), ithr
        real(dp),                 intent(in)    :: pens(:)
        real     :: lm(self%nrefs), w(1), neff, tins
        real(dp) :: lnz_k
        real     :: dmin_k
        integer  :: top(self%kin), inds(1), k, ri, nvalid, nent
        self%pen(:,i) = pens(1:self%nrefs)
        self%temp(i)  = posterior_temperature(pftc%wsqsum_ptcl_now(iptcl), pftc%get_pftsz())
        tins          = inpl_insert_temp(pftc, iptcl)
        nvalid = 0
        do ri = 1,self%nrefs
            if( self%dmin(ri,i) < 0.1 * huge(1.) )then
                lm(ri) = -marginal_dist(self%lnz(ri,i) - self%pen(ri,i), self%dmin(ri,i), self%temp(i))
                nvalid = nvalid + 1
            else
                lm(ri) = -huge(1.)
            endif
        enddo
        self%nin(i) = min(self%kin, nvalid)
        if( self%nin(i) < 1 ) return
        top = maxnloc(lm, self%kin)
        do k = 1,self%nin(i)
            ri = top(k)
            self%refs(k,i) = ri
            call inpl_row_at(self, pftc, iref_of(ri), iptcl, self%sh(:,ri,i), ithr)
            call inpl_quadrature(self%nrots, self%wrow(:,ithr), tins, 1, 0.d0, nent, inds, w, neff, lnz=lnz_k, &
                &dmin=dmin_k, post=self%rows(:,k,i))
        enddo
    end subroutine finish_row

    !> row i's marginals into marg (prob_row_marginal)
    subroutine get_marginal( self, i, marg )
        class(eul_prob_marginal), intent(in)    :: self
        integer,                  intent(in)    :: i
        type(prob_row_marginal),  intent(inout) :: marg
        if( .not. self%exists ) THROW_HARD('no marginal batch; eul_prob_marginal%get_marginal')
        if( i < self%i_first .or. i > self%i_last ) THROW_HARD('the row is not in the batch; eul_prob_marginal%get_marginal')
        call marg%kill
        marg%nrefs = self%nrefs
        marg%nrots = self%nrots
        marg%nin   = self%nin(i)
        marg%temp  = self%temp(i)
        allocate(marg%lnz(self%nrefs), source=self%lnz(:,i))
        allocate(marg%pen(self%nrefs), source=self%pen(:,i))
        allocate(marg%dmin(self%nrefs), source=self%dmin(:,i))
        allocate(marg%refs(marg%nin), source=self%refs(1:marg%nin,i))
        allocate(marg%rows(self%nrots, marg%nin), source=self%rows(:,1:marg%nin,i))
    end subroutine get_marginal

    !> The audit of entry ri of row i: the in-plane posterior insertion's own weighing of reference iref
    !! for particle iptcl at the shift sh_ins the insertion would take the entry at, with its shift-prior term pen_ins
    !! (the insertion's route: the particle shifted once above SHERRSQ, else the shifted-reference kernel; its
    !! temperature; inpl_quadrature at its budget), against the teacher's: dev(1:5) the absolute differences of the
    !! shift, the temperature, lnz, dmin and pen, dev(6) the largest of the in-plane posterior when the entry is one
    !! of the row's kept rows (0 otherwise), dev(7) 1 when it is kept; lnz_ins and dmin_ins the insertion's. On thread
    !! ithr's work
    subroutine audit_entry( self, pftc, iref, ri, i, iptcl, sh_ins, pen_ins, ithr, dev, lnz_ins, dmin_ins )
        class(eul_prob_marginal), intent(inout) :: self
        class(polarft_calc),      intent(inout) :: pftc
        integer,                  intent(in)    :: iref, ri, i, iptcl, ithr
        real,                     intent(in)    :: sh_ins(2)
        real(dp),                 intent(in)    :: pen_ins
        real(dp),                 intent(out)   :: dev(7), lnz_ins
        real,                     intent(out)   :: dmin_ins
        real     :: dists(self%nrots), post(self%nrots), w(INPL_POSTERIOR_BUDGET), neff, temp, dmin
        real(dp) :: lnz
        integer  :: inds(INPL_POSTERIOR_BUDGET), nent, k
        if( sum(sh_ins**2) > SHERRSQ )then
            call pftc%memoize_shifted_ptcl(iptcl, sh_ins, self%bank(:,:,ithr))
            call pftc%gen_euclid_dists_bank(iref, iptcl, self%bank(:,:,ithr), dists)
        else
            call pftc%gen_euclid_dists(iref, iptcl, sh_ins, dists)
        endif
        temp = temperature_or_argmin(posterior_temperature(pftc%wsqsum_ptcl_now(iptcl), pftc%get_pftsz()))
        call inpl_quadrature(self%nrots, dists, temp, INPL_POSTERIOR_BUDGET, 0.5d0, nent, inds, w, neff, lnz=lnz, &
            &dmin=dmin, post=post)
        if( nent < 1 ) dmin = huge(1.)
        dev    = 0.d0
        dev(1) = maxval(abs(real(sh_ins,dp) - real(self%sh(:,ri,i),dp)))
        dev(2) = abs(real(temp,dp) - real(temperature_or_argmin(self%temp(i)),dp))
        dev(3) = abs(lnz - self%lnz(ri,i))
        dev(4) = abs(real(dmin,dp) - real(self%dmin(ri,i),dp))
        dev(5) = abs(pen_ins - self%pen(ri,i))
        do k = 1,self%nin(i)
            if( self%refs(k,i) /= ri ) cycle
            dev(6) = maxval(abs(real(post,dp) - real(self%rows(:,k,i),dp)))
            dev(7) = 1.d0
        enddo
        lnz_ins  = lnz
        dmin_ins = dmin
    end subroutine audit_entry

    !> The audit of row i's class marginal: the insertion's view weights (inpl_view_weights, no extra
    !! dropped) over the row's kept references from the insertion's terms lnz_ins, dmin_ins and pen_ins (nrefs), against
    !! the teacher's row losses drow (the marginal's loss, as the rows file carries it) weighed at the physical
    !! temperature over the same references: dev(1) the largest absolute difference of the weights, dev(2) the bound
    !! of the losses' single-precision rounding on it
    subroutine audit_view_weights( self, i, lnz_ins, dmin_ins, pen_ins, drow, dev )
        class(eul_prob_marginal), intent(in)  :: self
        integer,                  intent(in)  :: i
        real(dp),                 intent(in)  :: lnz_ins(:), pen_ins(:)
        real,                     intent(in)  :: dmin_ins(:), drow(:)
        real(dp),                 intent(out) :: dev(2)
        real(dp) :: lm(self%kin), mt(self%kin)
        real     :: m(self%kin), temp
        integer  :: k, n
        dev = 0.d0
        n   = self%nin(i)
        if( n < 1 ) return
        temp = temperature_or_argmin(self%temp(i))
        call inpl_view_weights(n, lnz_ins(self%refs(1:n,i)), dmin_ins(self%refs(1:n,i)), temp, &
            &pen_ins(self%refs(1:n,i)), 0., m(1:n))
        do k = 1,n
            lm(k) = -real(temp,dp) * real(drow(self%refs(k,i)),dp)
        enddo
        mt(1:n) = exp(lm(1:n) - maxval(lm(1:n)))
        mt(1:n) = mt(1:n) / sum(mt(1:n))
        dev(1)  = maxval(abs(real(m(1:n),dp) - mt(1:n)))
        ! a loss carries a rounding of half its spacing; the weights move by at most temp times twice that
        dev(2)  = 2.d0 * real(temp,dp) * real(maxval(spacing(abs(drow(self%refs(1:n,i))))),dp) + 1.d-6
    end subroutine audit_view_weights

    !> reference iref's in-plane row against particle iptcl at shift into thread ithr's row, by the in-plane insertions'
    !! route (the particle shifted once above SHERRSQ, else the shifted-reference kernel)
    subroutine inpl_row_at( self, pftc, iref, iptcl, shift, ithr )
        class(eul_prob_marginal), intent(inout) :: self
        class(polarft_calc),      intent(inout) :: pftc
        integer,                  intent(in)    :: iref, iptcl, ithr
        real,                     intent(in)    :: shift(2)
        if( sum(shift**2) > SHERRSQ )then
            call pftc%memoize_shifted_ptcl(iptcl, shift, self%bank(:,:,ithr))
            call pftc%gen_euclid_dists_bank(iref, iptcl, self%bank(:,:,ithr), self%wrow(:,ithr))
        else
            call pftc%gen_euclid_dists(iref, iptcl, shift, self%wrow(:,ithr))
        endif
    end subroutine inpl_row_at

    !> thread ithr's row weighed for entry ri of row i at temp by the insertion's quadrature: lnz and dmin (the marker
    !! for no valid distance)
    subroutine weigh_row( self, ri, i, temp, ithr )
        class(eul_prob_marginal), intent(inout) :: self
        integer,                  intent(in)    :: ri, i, ithr
        real,                     intent(in)    :: temp
        real    :: w(1), neff
        integer :: inds(1), nent
        call inpl_quadrature(self%nrots, self%wrow(:,ithr), temp, 1, 0.d0, nent, inds, w, neff, lnz=self%lnz(ri,i), &
            &dmin=self%dmin(ri,i))
        if( nent < 1 ) self%dmin(ri,i) = huge(1.)
    end subroutine weigh_row

    subroutine kill( self )
        class(eul_prob_marginal), intent(inout) :: self
        if( allocated(self%lnz)  ) deallocate(self%lnz, self%pen, self%dmin, self%sh, self%temp, self%nin, self%refs, self%rows)
        if( allocated(self%wrow) ) deallocate(self%wrow, self%wcdf, self%bank)
        self%nrefs   = 0
        self%nrots   = 0
        self%kin     = 0
        self%i_first = 0
        self%i_last  = -1
        self%exists  = .false.
    end subroutine kill

end module simple_eul_prob_marginal
