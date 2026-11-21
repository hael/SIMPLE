!@descr: polarft class submodule: the block kernel, the posterior in-plane draws of a block of particles against many
!! references as one complex matrix product per angular bin, with its per-thread work and inverse transform plan
submodule (simple_polarft_calc) simple_polarft_block
use simple_linalg, only: cgemm
use simple_rnd,    only: rnd_keyed, rnd_draw_key, RND_FAM_INPL
#include "simple_local_flags.inc"
implicit none

complex(sp), parameter :: CONE = cmplx(1.,0.,kind=sp), CZERO = cmplx(0.,0.,kind=sp)

contains

    !> the work of a team of nthr_team and the plan of one particle's tile (PFTC_BLOCK_NR inverse transforms, out of
    !! place), planned with the calculator's other plans
    module subroutine new_block_work( self, nthr_team, plan_flags )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: nthr_team
        integer(c_int),      intent(in)    :: plan_flags
        integer :: ithr, n1, npairs
        call self%kill_block_work
        n1     = self%pftsz + 1
        npairs = PFTC_BLOCK_NR * PFTC_BLOCK_NB
        allocate(self%bwork(nthr_team))
        do ithr = 1,nthr_team
            ! not touched here (the kernel writes what it reads): a calculator that never scores a block touches only the
            ! pages its planner measures on
            allocate(self%bwork(ithr)%pterms(2*self%nk,PFTC_BLOCK_NB,n1), self%bwork(ithr)%tile(2*self%nk,PFTC_BLOCK_NR),&
                &self%bwork(ithr)%prods(PFTC_BLOCK_NB,PFTC_BLOCK_NR,n1), self%bwork(ithr)%norms(PFTC_BLOCK_NB))
            self%bwork(ithr)%pct  = fftwf_alloc_complex(int(n1 * npairs, c_size_t))
            call c_f_pointer(self%bwork(ithr)%pct, self%bwork(ithr)%ct, [n1, npairs])
            self%bwork(ithr)%pout = fftwf_alloc_real(int(self%nrots * npairs, c_size_t))
            call c_f_pointer(self%bwork(ithr)%pout, self%bwork(ithr)%out, [self%nrots, npairs])
        enddo
        self%plan_block = fftwf_plan_many_dft_c2r(1_c_int, [int(self%nrots,c_int)], int(PFTC_BLOCK_NR,c_int),&
            &self%bwork(1)%ct, [int(n1,c_int)], 1_c_int, int(n1,c_int), self%bwork(1)%out, [int(self%nrots,c_int)],&
            &1_c_int, int(self%nrots,c_int), plan_flags)
        if( .not. c_associated(self%plan_block) ) THROW_HARD('FFTW failed to create plan_block')
    end subroutine new_block_work

    module subroutine kill_block_work( self )
        class(polarft_calc), intent(inout) :: self
        integer :: ithr
        if( allocated(self%bwork) )then
            do ithr = 1,size(self%bwork)
                call fftwf_free(self%bwork(ithr)%pct)
                call fftwf_free(self%bwork(ithr)%pout)
            enddo
            deallocate(self%bwork)
        endif
        if( c_associated(self%plan_block) ) call fftwf_destroy_plan(self%plan_block)
        self%plan_block = c_null_ptr
    end subroutine kill_block_work

    !> gen_posterior_val at zero shift (gen_posterior_val_bank where l_bank) for nb particles of one half against nr
    !! references: per tile of PFTC_BLOCK_NR references both sums of every angular bin are one product over the 2nk
    !! shell terms (cgemm), then each pair's inverse transform and the per-pair kernels' distances and keyed draw. The
    !! product always has the shape PFTC_BLOCK_NB x PFTC_BLOCK_NR, a particle on any row (absent ones zero) and a
    !! reference on the column of its place in irefs, so that a pair's value depends on neither its block nor its
    !! thread, as the BLAS gives it. The sums (weights per term, their order) and the inverse transform (batched, out of
    !! place) differ from the per-pair kernels' in the last bits only. Thread-safe: the calling thread's work.
    !! objfun=euclid only.
    module subroutine score_block_posterior( self, nb, iptcls, nr, irefs, temps, dists, irots, banks, l_bank )
        class(polarft_calc),   intent(inout) :: self
        integer,               intent(in)    :: nb, nr
        integer,               intent(in)    :: iptcls(nb), irefs(nr)
        real(sp),              intent(in)    :: temps(nb)
        real(sp),              intent(out)   :: dists(nr,nb)
        integer,               intent(out)   :: irots(nr,nb)
        complex(sp), optional, intent(in)    :: banks(:,:,:)
        logical,     optional, intent(in)    :: l_bank(nb)
        real(dp) :: u
        real(sp) :: wk
        integer  :: ithr, i, i1, ih, ieo, b, j, k, kk, k0, n1, n2k, r0, nt, ir, iref0, icol
        logical  :: l_banks, l_consec, l_b
        call self%check_ptcl_memo
        call self%check_euclid_weights
        if( self%p_ptr%cc_objfun /= OBJFUN_EUCLID ) THROW_HARD('score_block_posterior requires objfun=euclid')
        if( .not. self%l_rnd_scope ) THROW_HARD('the posterior draw needs the calculator''s scope; score_block_posterior')
        if( self%ref_gen_bins /= self%ref_gen ) THROW_HARD('stale or missing bin-major memo: memoize_ref_bins after memoize_refs')
        if( nb < 1 .or. nb > PFTC_BLOCK_NB ) THROW_HARD('a block holds 1 to PFTC_BLOCK_NB particles; score_block_posterior')
        if( nr < 1 ) return
        if( minval(irefs) < 1 .or. maxval(irefs) > self%nrefs ) THROW_HARD('reference out of range; score_block_posterior')
        if( present(banks) .neqv. present(l_bank) ) THROW_HARD('banks and l_bank go together; score_block_posterior')
        l_banks = present(banks)
        if( l_banks )then
            if( size(banks,1) /= self%pftsz+1 .or. size(banks,2) /= self%nk .or. size(banks,3) < nb )then
                THROW_HARD('the banks have the wrong shape; score_block_posterior')
            endif
        endif
        ithr = omp_get_thread_num() + 1
        if( ithr > size(self%bwork) ) THROW_HARD('the team exceeds the block kernel''s work; score_block_posterior')
        n1  = self%pftsz + 1
        n2k = 2 * self%nk
        k0  = self%kfromto(1)
        i1  = self%pinds(iptcls(1))
        ih  = merge(1, 2, self%iseven(i1))
        ieo = merge(REF_EVEN, REF_ODD, self%iseven(i1))
        associate( bw => self%bwork(ithr) )
            ! the particle terms of every bin: rows 1..nk w_k FT(CTF2), rows nk+1..2nk -2 w_k FT(X.CTF), zero rows past nb
            do b = 1,nb
                i = self%pinds(iptcls(b))
                if( self%iseven(i) .neqv. self%iseven(i1) ) THROW_HARD('a block mixes the halves; score_block_posterior')
                bw%norms(b) = self%wsqsum_ptcl_now(iptcls(b)) * real(2*self%nrots, dp)
                l_b = .false.
                if( l_banks ) l_b = l_bank(b)
                do j = 1,n1
                    do k = self%kfromto(1),self%kfromto(2)
                        kk = k - k0 + 1
                        wk = self%sigma2_wk(k,ih)
                        bw%pterms(kk,b,j) = wk * self%ft_ctf2(j,k,i)
                    enddo
                enddo
                if( l_b )then
                    do j = 1,n1
                        do kk = 1,self%nk
                            bw%pterms(self%nk+kk,b,j) = (-2._sp * self%sigma2_wk(k0+kk-1,ih)) * banks(j,kk,b)
                        enddo
                    enddo
                else
                    do j = 1,n1
                        do k = self%kfromto(1),self%kfromto(2)
                            kk = k - k0 + 1
                            bw%pterms(self%nk+kk,b,j) = (-2._sp * self%sigma2_wk(k,ih)) * self%ft_ptcl_ctf(j,k,i)
                        enddo
                    enddo
                endif
            enddo
            if( nb < PFTC_BLOCK_NB ) bw%pterms(:,nb+1:PFTC_BLOCK_NB,:) = CZERO
            do r0 = 1,nr,PFTC_BLOCK_NR
                nt    = min(PFTC_BLOCK_NR, nr - r0 + 1)
                iref0 = irefs(r0)
                l_consec = nt == PFTC_BLOCK_NR
                do ir = 2,nt
                    if( irefs(r0+ir-1) /= iref0 + ir - 1 ) l_consec = .false.
                enddo
                ! D(b,ir) = sum_k conj(P(k,b)) G(k,ir), the conjugate of the pair's sums: a full tile of consecutive
                ! references is read in place, any other tile gathered bin by bin with zero columns past nt
                if( l_consec )then
                    do j = 1,n1
                        call cgemm('C', 'N', PFTC_BLOCK_NB, PFTC_BLOCK_NR, n2k, CONE, bw%pterms(1,1,j), n2k,&
                            &self%ft_ref_bins(1,iref0,j,ieo), n2k, CZERO, bw%prods(1,1,j), PFTC_BLOCK_NB)
                    enddo
                else
                    if( nt < PFTC_BLOCK_NR ) bw%tile(:,nt+1:PFTC_BLOCK_NR) = CZERO
                    do j = 1,n1
                        do ir = 1,nt
                            bw%tile(:,ir) = self%ft_ref_bins(:,irefs(r0+ir-1),j,ieo)
                        enddo
                        call cgemm('C', 'N', PFTC_BLOCK_NB, PFTC_BLOCK_NR, n2k, CONE, bw%pterms(1,1,j), n2k,&
                            &bw%tile, n2k, CZERO, bw%prods(1,1,j), PFTC_BLOCK_NB)
                    enddo
                endif
                ! each pair's bins contiguous (the conjugate back), one particle's tile per inverse transform
                do b = 1,nb
                    do ir = 1,PFTC_BLOCK_NR
                        icol = (b-1)*PFTC_BLOCK_NR + ir
                        bw%ct(1:n1,icol) = conjg(bw%prods(b,ir,1:n1))
                    enddo
                    icol = (b-1)*PFTC_BLOCK_NR + 1
                    call fftwf_execute_dft_c2r(self%plan_block, bw%ct(:,icol:icol+PFTC_BLOCK_NR-1),&
                        &bw%out(:,icol:icol+PFTC_BLOCK_NR-1))
                enddo
                do b = 1,nb
                    do ir = 1,nt
                        icol = (b-1)*PFTC_BLOCK_NR + ir
                        u = rnd_keyed(rnd_draw_key(self%rnd_scope, iptcls(b), RND_FAM_INPL, reference=irefs(r0+ir-1)))
                        call draw_posterior_rvec(self, ithr, bw%out(:,icol), bw%norms(b), temps(b), u,&
                            &dists(r0+ir-1,b), irots(r0+ir-1,b))
                    enddo
                enddo
            enddo
        end associate
    end subroutine score_block_posterior

end submodule simple_polarft_block
