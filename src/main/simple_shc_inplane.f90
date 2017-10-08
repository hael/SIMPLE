module simple_shc_inplane
use simple_defs      ! use all in there
use simple_ran_tabu, only: ran_tabu
implicit none

type shc_inplane
    private
    type(ran_tabu)       :: rt               !< random order generator
    real,    allocatable :: shifts(:,:,:)    !< shift grid (NRADIAL_LINES,NK,2)
    logical              :: exists = .false. !< to flag existence
contains 
    procedure :: new
    procedure :: srch
    procedure :: kill
end type shc_inplane

integer, parameter :: NRADIAL_LINES = 32, KVECSZ = 8
real,    parameter :: KVEC(KVECSZ) = [0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0]

contains

    subroutine new( self )
        use simple_math,   only: gen_polar_coords
        class(shc_inplane), intent(inout) :: self
        real, allocatable :: angtab(:)
        call self%kill
        ! rotational origin shifts are searched on a polar grid, since it 
        ! biases the distribution toward the center, the rings closest to
        ! the center are evaluated first
        call gen_polar_coords(KVEC, NRADIAL_LINES, self%shifts, angtab)
        ! make random number generator
        self%rt = ran_tabu(NRADIAL_LINES)
        ! flag existence
        self%exists = .true.
    end subroutine new

    subroutine srch( self, pftcc, ref, iptcl, nrots_tot, prev_rot, rot, shvec )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(shc_inplane),      intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: ref, iptcl, nrots_tot, prev_rot
        integer,                 intent(out)   :: rot
        real,                    intent(out)   :: shvec(2)
        integer :: srch_order(NRADIAL_LINES), k, j, jj, loc(1)
        real    :: cc_prev, cc, corrs(nrots_tot)
        ! init
        rot     = prev_rot
        shvec   = [0.,0.]
        corrs   = pftcc%gencorrs_fft(ref, iptcl, shvec)
        cc_prev = maxval(corrs)
        ! srch 
        do k=1,KVECSZ
            call self%rt%reset
            call self%rt%ne_ran_iarr(srch_order)
            do j=1,NRADIAL_LINES
                jj    = srch_order(j)
                corrs = pftcc%gencorrs_fft(ref, iptcl, self%shifts(jj,k,:))
                loc   = maxloc(corrs)
                cc    = corrs(loc(1))
                if( cc > cc_prev )then
                    rot     = loc(1)
                    shvec   = self%shifts(jj,k,:)
                    cc_prev = cc
                    if( k >= 2 ) return ! minimum of half a pixels shift evaluated
                endif
            end do
        end do
    end subroutine srch

    subroutine kill( self )
        class(shc_inplane), intent(inout) :: self
        if( self%exists )then
            call self%rt%kill
            deallocate(self%shifts)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_shc_inplane
