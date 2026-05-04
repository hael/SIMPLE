!@descr: polarft class submodule for geometry-related things: shift, rotate, mirror etc.
submodule (simple_polarft_calc) simple_polarft_geom
implicit none

contains

    !>  Generate polar shift matrix by means of de Moivre's formula
    module subroutine gen_shmat(self, ithr, shift, shmat)
        class(polarft_calc),  intent(inout) :: self
        integer,              intent(in)    :: ithr
        real(sp),             intent(in)    :: shift(2)
        complex(sp), pointer, intent(inout) :: shmat(:,:)
        call self%gen_shmat_8(ithr, real(shift,dp), self%heap_vars(ithr)%shmat_8)
        shmat = cmplx(self%heap_vars(ithr)%shmat_8)
    end subroutine gen_shmat

    !>  Generate polar shift matrix by means of de Moivre's formula, double precision
    module subroutine gen_shmat_8(self, ithr, shift_8, shmat_8)
        class(polarft_calc),  intent(inout) :: self
        integer,              intent(in)    :: ithr
        real(dp),             intent(in)    :: shift_8(2)
        complex(dp), pointer, intent(inout) :: shmat_8(:,:)
         integer :: k
        ! first shell, analytic
        self%heap_vars(ithr)%argvec = self%argtransf(:self%pftsz,  self%kfromto(1)) * shift_8(1) +&
                                    & self%argtransf(self%pftsz+1:,self%kfromto(1)) * shift_8(2)
        shmat_8(:,self%kfromto(1))  = dcmplx(dcos(self%heap_vars(ithr)%argvec), dsin(self%heap_vars(ithr)%argvec))
        ! one shell to the next
        self%heap_vars(ithr)%argvec = self%argtransf_shellone(:self%pftsz)   * shift_8(1) +&
                                    & self%argtransf_shellone(self%pftsz+1:) * shift_8(2)
        self%heap_vars(ithr)%shvec  = dcmplx(dcos(self%heap_vars(ithr)%argvec), dsin(self%heap_vars(ithr)%argvec))
        ! remaining shells, cos(kx)+isin(kx) = (cos(x)+isin(x))**k-1 * (cos(x)+isin(x))
        do k = self%kfromto(1)+1,self%interpklim
            shmat_8(:,k) = shmat_8(:,k-1) * self%heap_vars(ithr)%shvec
        enddo
        ! alternative to:
        ! argmat  => self%heap_vars(ithr)%argmat_8
        ! argmat  =  self%argtransf(:self%pftsz,:)*shvec(1) + self%argtransf(self%pftsz + 1:,:)*shvec(2)
        ! shmat   =  cmplx(cos(argmat),sin(argmat),dp)
    end subroutine gen_shmat_8

    module subroutine shift_ptcl( self, iptcl, shvec)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        real(sp),            intent(in)    :: shvec(2)
        complex(sp), pointer :: shmat(:,:)
        integer :: ithr, i
        ithr  = omp_get_thread_num() + 1
        i     = self%pinds(iptcl)
        shmat => self%heap_vars(ithr)%shmat
        call self%gen_shmat(ithr, shvec, shmat)
        self%pfts_ptcls(:,:,i) = self%pfts_ptcls(:,:,i) * shmat
    end subroutine shift_ptcl

    module subroutine shift_ref( self, iref, shvec)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref
        real(sp),            intent(in)    :: shvec(2)
        complex(dp), pointer :: shmat(:,:)
        integer :: ithr
        ithr = omp_get_thread_num() + 1
        shmat => self%heap_vars(ithr)%shmat_8
        call self%gen_shmat_8(ithr, real(shvec,dp), shmat)
        self%pfts_refs_even(:,:,iref) = cmplx(dcmplx(self%pfts_refs_even(:,:,iref)) * shmat)
        self%pfts_refs_odd( :,:,iref) = cmplx(dcmplx(self%pfts_refs_odd( :,:,iref)) * shmat)
    end subroutine shift_ref

    ! mirror pft about h (mirror about y of cartesian image)
    module pure subroutine mirror_ref_pft( self, pft, pftm )
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
    end subroutine mirror_ref_pft

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
