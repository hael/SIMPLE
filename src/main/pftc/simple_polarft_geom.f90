!@descr: polarft class submodule for geometry-related things: shift, rotate, mirror etc.
submodule (simple_polarft_calc) simple_polarft_geom
#include "simple_local_flags.inc"
implicit none

contains

    !> Mirror one bank into its canonical frames; general angles use trigonometric interpolation.
    module subroutine mirror_ref_bank(self, refs, partners, frame_angles)
        use simple_ori_utils, only: canonical_mirror_frame_angle
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        class(polarft_calc), intent(inout) :: self
        complex(sp),        intent(inout) :: refs(:,:,:)
        integer,            intent(in)    :: partners(:)
        real,               intent(in)    :: frame_angles(:)
        integer :: iproj, m, ir, jr, kk, ithr, mode, freq, nhalf
        real(dp) :: psi, arg
        complex(sp) :: phase
        logical :: half_turn
        nhalf = size(partners)
        if( OMP_IN_PARALLEL() ) THROW_HARD('No bank mirroring inside OpenMP regions')
        if( .not. self%existence ) THROW_HARD('bank mirroring requires an initialized polar calculator')
        if( size(refs,1) /= self%pftsz .or. size(refs,2) /= self%nk ) THROW_HARD('mirror bank polar dimensions differ')
        if( size(refs,3) /= 2*nhalf .or. size(frame_angles) /= nhalf ) THROW_HARD('mirror bank pair dimensions differ')
        if( .not. all(ieee_is_finite(frame_angles)) ) THROW_HARD('mirror frame angles must be finite')
        if( omp_get_max_threads() > size(self%cmat2_many) ) THROW_HARD('mirror bank team exceeds its workspace')
        !$omp parallel do default(shared) private(iproj,m,ir,jr,kk,ithr,mode,freq,psi,arg,phase,half_turn) schedule(static) proc_bind(close)
        do iproj = 1,nhalf
            ithr = omp_get_thread_num()+1
            m = partners(iproj)
            psi = canonical_mirror_frame_angle(frame_angles(iproj))
            half_turn = psi == 180._dp
            if( psi == 0._dp .or. half_turn )then
                do kk = 1,self%nk
                    refs(1,kk,m) = merge(refs(1,kk,iproj),conjg(refs(1,kk,iproj)),half_turn)
                    do ir = 2,self%pftsz
                        jr = self%pftsz-ir+2
                        refs(ir,kk,m) = merge(conjg(refs(jr,kk,iproj)),refs(jr,kk,iproj),half_turn)
                    enddo
                enddo
            else
                self%cmat2_many(ithr)%c(1:self%pftsz,:) = refs(:,:,iproj)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,:) = conjg(refs(:,:,iproj))
                call fftwf_execute_dft(self%plan_fwd1_many,self%cmat2_many(ithr)%c,self%cmat2_many(ithr)%c)
                psi = psi*DPI/180._dp
                do mode = 0,self%nrots-1
                    freq = mode
                    if( mode > self%pftsz ) freq = mode-self%nrots
                    arg = real(freq,dp)*psi
                    phase = cmplx(cos(arg),sin(arg),sp)
                    ! Split the aliased +/- Nyquist modes equally to preserve Hermitian half-turn symmetry.
                    if( mode == self%pftsz ) phase = cmplx(cos(arg),0._dp,sp)
                    self%cmat2_many(ithr)%c(mode+1,:) = self%cmat2_many(ithr)%c(mode+1,:)*phase
                enddo
                ! conjg(FFT(A*exp(+im*psi)))/n evaluates the reflected ring in the zero-e3 frame.
                call fftwf_execute_dft(self%plan_fwd1_many,self%cmat2_many(ithr)%c,self%cmat2_many(ithr)%c)
                refs(:,:,m) = conjg(self%cmat2_many(ithr)%c(1:self%pftsz,:))/real(self%nrots,sp)
            endif
        enddo
        !$omp end parallel do
    end subroutine mirror_ref_bank

    !>  Generate polar shift matrix by means of de Moivre's formula, for alignment: only the searched shells
    !>  [kfromto(1);kfromto(2)] are set, the only ones its consumers read
    module subroutine gen_shmat4aln(self, ithr, shift, shmat)
        class(polarft_calc),  intent(inout) :: self
        integer,              intent(in)    :: ithr
        real(sp),             intent(in)    :: shift(2)
        complex(sp), pointer, intent(inout) :: shmat(:,:)
        call self%gen_shmat4aln_8(ithr, real(shift,dp), self%heap_vars(ithr)%shmat_8)
        shmat(:,self%kfromto(1):self%kfromto(2)) = cmplx(self%heap_vars(ithr)%shmat_8(:,self%kfromto(1):self%kfromto(2)))
    end subroutine gen_shmat4aln

    !>  Generate polar shift matrix by means of de Moivre's formula, double precision
    !>  To be used only for alignment as it is band-passed within [kfromto(1);kfromto(2)]
    module subroutine gen_shmat4aln_8(self, ithr, shift_8, shmat_8)
        class(polarft_calc),  intent(inout) :: self
        integer,              intent(in)    :: ithr
        real(dp),             intent(in)    :: shift_8(2)
        complex(dp), pointer, intent(inout) :: shmat_8(:,:)
        real(dp) :: rarg
        integer  :: k,p
        ! per shell multiplicative increment
        do p = 1,self%pftsz
            rarg =        self%argtransf_shellone(           p) * shift_8(1)
            rarg = rarg + self%argtransf_shellone(self%pftsz+p) * shift_8(2)
            self%heap_vars(ithr)%shvec(p) = cmplx(cos(rarg), sin(rarg), kind=dp)
        enddo
        ! first shell
        k = self%kfromto(1)
        do p = 1,self%pftsz
            rarg         =        self%argtransf(           p,k) * shift_8(1)
            rarg         = rarg + self%argtransf(self%pftsz+p,k) * shift_8(2)
            shmat_8(p,k) = cmplx(cos(rarg), sin(rarg), kind=dp)
        enddo
        ! remaining shells, cos(kx)+isin(kx) = (cos(x)+isin(x))**k-1 * (cos(x)+isin(x))
        do k = self%kfromto(1)+1,self%kfromto(2)
            do p = 1,self%pftsz
                shmat_8(p,k) = shmat_8(p,k-1) * self%heap_vars(ithr)%shvec(p)
            enddo
        enddo
    end subroutine gen_shmat4aln_8

    ! mirror pft about h (mirror about y of cartesian image)
    module pure subroutine mirror_ref_pft_by_pft( self, pft, pftm )
        class(polarft_calc), intent(in)  :: self
        complex(sp),         intent(in)  :: pft(1:self%pftsz, self%kfromto(1):self%interpklim)
        complex(sp),         intent(out) :: pftm(1:self%pftsz, self%kfromto(1):self%interpklim)
        integer :: i,j
        pftm(1,:) = conjg(pft(1,:))
        if( is_even(self%pftsz) )then
            do i = 2,self%pftsz/2
                j = self%pftsz-i+2
                pftm(i,:) = pft(j,:)
                pftm(j,:) = pft(i,:)
            enddo
            i = self%pftsz/2 + 1
            pftm(i,:) = pft(i,:)
        else
            do i = 2,(self%pftsz+1)/2
                j = self%pftsz-i+2
                pftm(i,:) = pft(j,:)
                pftm(j,:) = pft(i,:)
            enddo
        endif
    end subroutine mirror_ref_pft_by_pft

    ! mirrors the even reference into the odd position
    module pure subroutine mirror_ref_pft_by_ind( self, iref )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref
        call self%mirror_ref_pft_by_pft(self%pfts_refs(:,:,iref,REF_EVEN), self%pfts_refs(:,:,iref,REF_ODD))
    end subroutine mirror_ref_pft_by_ind

    module subroutine rotate_ref_8( self, pft, irot, pft_rot)
        class(polarft_calc), intent(in)  :: self
        complex(dp),         intent(in)  :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,             intent(in)  :: irot
        complex(dp),         intent(out) :: pft_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer :: mid
        if( irot == 1 )then
            pft_rot = pft
        elseif( irot >= 2 .and. irot <= self%pftsz )then
            mid = self%pftsz - irot + 1
            pft_rot(   1:irot-1,    :) = conjg(pft(mid+1:self%pftsz,:))
            pft_rot(irot:self%pftsz,:) =       pft(    1:mid,       :)
        elseif( irot == self%pftsz + 1 )then
            pft_rot = conjg(pft)
        else
            mid = self%nrots - irot + 1
            pft_rot(irot-self%pftsz:self%pftsz       ,:) = conjg(pft(1    :mid,       :))
            pft_rot(1              :irot-self%pftsz-1,:) =       pft(mid+1:self%pftsz,:)
        endif
    end subroutine rotate_ref_8

    ! private particle rotation utility
    module subroutine rotate_ptcl(self, i, irot, pft_rot)
        class(polarft_calc), intent(in)  :: self
        integer,             intent(in)  :: i       ! internal particle index
        integer,             intent(in)  :: irot
        complex(sp),         intent(out) :: pft_rot(self%pftsz,self%kfromto(1):self%interpklim)
        integer :: mid
        if( irot == 1 )then
            pft_rot = self%pfts_ptcls(:,:,i)
        elseif( irot >= 2 .and. irot <= self%pftsz )then
            mid = self%pftsz - irot + 1
            pft_rot(   1:irot-1,    :) = conjg(self%pfts_ptcls(mid+1:self%pftsz,:,i))
            pft_rot(irot:self%pftsz,:) =       self%pfts_ptcls(    1:mid,       :,i)
        elseif( irot == self%pftsz + 1 )then
            pft_rot = conjg(self%pfts_ptcls(:,:,i))
        else
            mid = self%nrots - irot + 1
            pft_rot(irot-self%pftsz:self%pftsz       ,:) = conjg(self%pfts_ptcls(1    :mid,       :,i))
            pft_rot(1              :irot-self%pftsz-1,:) =       self%pfts_ptcls(mid+1:self%pftsz,:,i)
        endif
    end subroutine rotate_ptcl

    ! private |CTF| rotation utility
    module subroutine rotate_ctf(self, i, irot, ctf_rot)
        class(polarft_calc), intent(in)  :: self
        integer,             intent(in)  :: i       ! internal particle index
        integer,             intent(in)  :: irot
        real(sp),            intent(out) :: ctf_rot(self%pftsz,self%kfromto(1):self%interpklim)
        integer :: mid
        if( irot == 1 )then
            ctf_rot = self%ctfmats(:,:,i)
        elseif( irot >= 2 .and. irot <= self%pftsz )then
            mid = self%pftsz - irot + 1
            ctf_rot(   1:irot-1,    :) = self%ctfmats(mid+1:self%pftsz,:,i)
            ctf_rot(irot:self%pftsz,:) = self%ctfmats(    1:mid,       :,i)
        elseif( irot == self%pftsz + 1 )then
            ctf_rot = self%ctfmats(:,:,i)
        else
            mid = self%nrots - irot + 1
            ctf_rot(irot-self%pftsz:self%pftsz       ,:) = self%ctfmats(1    :mid,       :,i)
            ctf_rot(1              :irot-self%pftsz-1,:) = self%ctfmats(mid+1:self%pftsz,:,i)
        endif
    end subroutine rotate_ctf

end submodule simple_polarft_geom
