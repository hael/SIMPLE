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

    ! Calculates weights for two neighboring lines based on linear Euclidean
    ! distance in the 2D Fourier plane. The weights are normalized so that lw + rw = 1 for each k.
    module pure subroutine gen_clin_weights( self, psi, lrot, rrot, lw, rw )
        class(polarft_calc), intent(in)    :: self
        real(dp),            intent(in)    :: psi
        integer,             intent(inout) :: lrot, rrot
        real(dp),            intent(out)   :: lw(self%kfromto(1):self%interpklim)
        real(dp),            intent(out)   :: rw(self%kfromto(1):self%interpklim)
        ! Local variables
        real(dp) :: sinpsi, cospsi       ! sin and cos of the target angle
        real(dp) :: h_line, k_line       ! Cartesian coords of the point on the target common line
        real(dp) :: l_h_polar, l_k_polar ! Cartesian coords of the point on the left reference line
        real(dp) :: r_h_polar, r_k_polar ! Cartesian coords of the point on the right reference line
        real(dp) :: dist_l, dist_r       ! Euclidean distance from target point to left/right points
        real(dp) :: dist_tot             ! Sum of the two distances
        real(dp) :: ang_l, d, eps        ! Temporary variables for angle calculations and numerical stability
        real(dp) :: kw                   ! Jacobian weight to account for the increasing density of points in Fourier space with increasing radius
        integer  :: k                    ! Resolution index
        ! --- 1. Find the neighboring reference angles ---
        lrot = self%get_roind_fast(real(psi))
        ang_l = real(self%angtab(lrot), dp)
        ! Simple difference check to find the correct neighbors
        d = psi - ang_l
        if (d >  180.0_dp) d = d - 360.0_dp
        if (d < -180.0_dp) d = d + 360.0_dp
        if (d < 0.0_dp) then
            ! psi is to the "left" of the closest angle, so the closest is the right neighbor
            rrot = lrot
            lrot = lrot - 1
            if (lrot < 1) lrot = lrot + self%get_nrots()
        else
            ! psi is to the "right", so the closest is the left neighbor
            rrot = lrot + 1
            if (rrot > self%get_nrots()) rrot = rrot - self%get_nrots()
        end if
        ! Pre-calculate sin and cos for the main loop
        sinpsi = sin(deg2rad(psi))
        cospsi = cos(deg2rad(psi))
        ! --- 2. Loop over all radial points (k) to calculate weights ---
        do k = self%kfromto(1), self%interpklim
            ! LEFT reference line point (h,k) in Cartesian coords
            l_h_polar = real(self%polar(1,k,lrot), dp)
            l_k_polar = real(self%polar(2,k,lrot), dp)
            ! RIGHT reference line point (h,k) in Cartesian coords
            r_h_polar = real(self%polar(1,k,rrot), dp)
            r_k_polar = real(self%polar(2,k,rrot), dp)
            !  TARGET line point (h,k) in Cartesian coords
            h_line    =  sinpsi * real(k, dp)
            k_line    = -cospsi * real(k, dp)
            ! Distances in the Fourier plane from the target point to the left and right reference points
            dist_l    = sqrt((l_h_polar - h_line)**2 + (l_k_polar - k_line)**2)
            dist_r    = sqrt((r_h_polar - h_line)**2 + (r_k_polar - k_line)**2)
            dist_tot  = dist_l + dist_r
            ! Robust "on-line" eps (scale-aware)
            eps       = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, dist_tot)
            if (dist_l <= eps .and. dist_r <= eps) then
                lw(k) = 0.5_dp
                rw(k) = 0.5_dp
            else if (dist_l <= eps) then
                lw(k) = 1.0_dp
                rw(k) = 0.0_dp
            else if (dist_r <= eps) then
                lw(k) = 0.0_dp
                rw(k) = 1.0_dp
            else
                lw(k) = dist_r / dist_tot
                rw(k) = dist_l / dist_tot
            end if
            kw = real(k, dp)
            lw(k) = lw(k) * kw
            rw(k) = rw(k) * kw
        end do
    end subroutine gen_clin_weights

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
