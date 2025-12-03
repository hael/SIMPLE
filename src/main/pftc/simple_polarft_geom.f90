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
        do k = self%kfromto(1)+1,self%kfromto(2)
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
    module subroutine mirror_ref_pft( self, iref )
        class(polarft_calc), target, intent(in) :: self
        integer,                     intent(in) :: iref
        integer :: i,j
        complex(sp), pointer :: pft(:,:) => null(), pftmirr(:,:) => null()
        pft     => self%pfts_refs_even(:,:,iref)
        pftmirr => self%pfts_refs_odd(:,:,iref)
        pftmirr(1,:) = conjg(pft(1,:))
        if( is_even(self%pftsz) )then
            do i = 2,self%pftsz/2
                j = self%pftsz-i+2
                pftmirr(i,:) = pft(j,:)
                pftmirr(j,:) = pft(i,:)
            enddo
            i = self%pftsz/2 + 1
            pftmirr(i,:) = pft(i,:)
        else
            do i = 2,(self%pftsz+1)/2
                j = self%pftsz-i+2
                pftmirr(i,:) = pft(j,:)
                pftmirr(j,:) = pft(i,:)
            enddo
        endif
    end subroutine mirror_ref_pft

    module subroutine rotate_pft_1(self, pft, irot, pft_rot)
        class(polarft_calc), intent(in)  :: self
        complex(dp),         intent(in)  :: pft(:,:)
        integer,             intent(in)  :: irot
        complex(dp),         intent(out) :: pft_rot(:,:)
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
    end subroutine rotate_pft_1

    module subroutine rotate_pft_2(self, pft, irot, pft_rot)
        class(polarft_calc), intent(in)  :: self
        complex(sp),         intent(in)  :: pft(:,:)
        integer,             intent(in)  :: irot
        complex(sp),         intent(out) :: pft_rot(:,:)
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
    end subroutine rotate_pft_2

    module subroutine rotate_pft_3(self, pft, irot, pft_rot)
        class(polarft_calc), intent(in)  :: self
        real(sp),            intent(in)  :: pft(:,:)
        integer,             intent(in)  :: irot
        real(sp),            intent(out) :: pft_rot(:,:)
        integer :: mid
        if( irot == 1 )then
            pft_rot = pft
        elseif( irot >= 2 .and. irot <= self%pftsz )then
            mid = self%pftsz - irot + 1
            pft_rot(   1:irot-1,    :) = pft(mid+1:self%pftsz,:)
            pft_rot(irot:self%pftsz,:) = pft(    1:mid,       :)
        elseif( irot == self%pftsz + 1 )then
            pft_rot = pft
        else
            mid = self%nrots - irot + 1
            pft_rot(irot-self%pftsz:self%pftsz       ,:) = pft(1    :mid,       :)
            pft_rot(1              :irot-self%pftsz-1,:) = pft(mid+1:self%pftsz,:)
        endif
    end subroutine rotate_pft_3

    module subroutine rotate_pft_4(self, pft, irot, pft_rot)
        class(polarft_calc), intent(in)  :: self
        real(dp),            intent(in)  :: pft(:,:)
        integer,             intent(in)  :: irot
        real(dp),            intent(out) :: pft_rot(:,:)
        integer :: mid
        if( irot == 1 )then
            pft_rot = pft
        elseif( irot >= 2 .and. irot <= self%pftsz )then
            mid = self%pftsz - irot + 1
            pft_rot(   1:irot-1,    :) = pft(mid+1:self%pftsz,:)
            pft_rot(irot:self%pftsz,:) = pft(    1:mid,       :)
        elseif( irot == self%pftsz + 1 )then
            pft_rot = pft
        else
            mid = self%nrots - irot + 1
            pft_rot(irot-self%pftsz:self%pftsz       ,:) = pft(1    :mid,       :)
            pft_rot(1              :irot-self%pftsz-1,:) = pft(mid+1:self%pftsz,:)
        endif
    end subroutine rotate_pft_4

    module subroutine rotate_iref_1(self, iref, irot, sh)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, irot
        real,                intent(in)    :: sh(2)
        complex(dp), pointer :: pft_ref(:,:), pft_ref_tmp_8(:,:), shmat(:,:)
        integer :: ithr
        ithr = omp_get_thread_num() + 1
        pft_ref       => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat         => self%heap_vars(ithr)%shmat_8
        ! even
        pft_ref_tmp_8 = dcmplx(self%pfts_refs_even(:,:,iref))
        call self%gen_shmat_8(ithr, real(sh,dp), shmat)
        call self%rotate_pft(pft_ref_tmp_8, irot, pft_ref)
        pft_ref                       = pft_ref * shmat
        self%pfts_refs_even(:,:,iref) = cmplx(pft_ref)
        ! odd
        pft_ref_tmp_8                 = dcmplx(self%pfts_refs_odd(:,:,iref))
        call self%rotate_pft(pft_ref_tmp_8, irot, pft_ref)
        pft_ref                       = pft_ref * shmat
        self%pfts_refs_odd( :,:,iref) = cmplx(pft_ref)
    end subroutine rotate_iref_1

    module subroutine rotate_iref_2(self, pft_ref, irot, sh, pft_ref_out)
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(in)    :: pft_ref(:,:)
        integer,              intent(in)    :: irot
        real,                 intent(in)    :: sh(2)
        complex(dp), pointer, intent(out)   :: pft_ref_out(:,:)
        complex(dp), pointer :: shmat(:,:)
        integer :: ithr
        ithr  = omp_get_thread_num() + 1
        shmat => self%heap_vars(ithr)%shmat_8
        call self%rotate_pft(pft_ref, irot, pft_ref_out)
        call self%gen_shmat_8(ithr, real(sh,dp), shmat)
        pft_ref_out = pft_ref_out * shmat
    end subroutine rotate_iref_2

end submodule simple_polarft_geom
