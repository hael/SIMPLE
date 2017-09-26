module simple_shc_inplane
use simple_defs      ! use all in there
use simple_ran_tabu, only: ran_tabu
implicit none

type shc_inplane
    private
    type(ran_tabu)       :: rt              !< random order generator
    integer, allocatable :: srch_order(:)   !< random search order
    integer, allocatable :: roinds(:)       !< rotation index grid
    real,    allocatable :: shifts(:,:)     !< shift grid
    integer              :: maxevals        !< max # cc evals
    integer              :: minevals        !< min # cc evals
    logical              :: exists=.false.  !< to flag existence
contains 
    procedure          :: new
    procedure, private :: update_grid
    procedure          :: srch
    procedure          :: kill
end type shc_inplane

contains

    subroutine new( self )
        use simple_syslib, only: alloc_errchk
        class(shc_inplane), intent(inout) :: self
        type(ran_tabu) :: rt
        integer        :: nrots_here, nsh, nall, alloc_stat
        call self%kill
        ! # rotations
        nrots_here = SHC_INPL_INPLHWDTH * 2 + 1
        ! # shifts
        nsh  = nint((SHC_INPL_TRSHWDTH * 2.0) / SHC_INPL_TRSSTEPSZ)
        nsh  = nsh * nsh
        ! # all, the -nrots_here comes from the prev best (prev_rot,0,0) config 
        nall = nrots_here * nsh - nrots_here
        ! maximum # cc evals is half of the search space
        self%maxevals = nall / 2
        ! minimum # cc evals is 10% of the search space
        self%minevals = nint(real(nall)*0.1)
        ! allocate
        allocate(self%srch_order(nall), self%roinds(nall), self%shifts(nall,2), stat=alloc_stat)
        call alloc_errchk('In: new; simple_shc_inplane', alloc_stat)
        self%srch_order = 0
        self%roinds     = 0
        self%shifts     = 0.0
        ! make random number generator
        self%rt = ran_tabu(nall)
        ! flag existence
        self%exists = .true.
    end subroutine new

    subroutine update_grid( self, nrots, prev_rot )
        use simple_math, only: cyci_1d
        class(shc_inplane), intent(inout) :: self
        integer,            intent(in)    :: nrots, prev_rot
        integer        :: j, jrot, cnt
        real           :: xsh, ysh
        ! SETUP THE OPTIMISATION GRID
        cnt = 0
        do j=prev_rot - SHC_INPL_INPLHWDTH,prev_rot + SHC_INPL_INPLHWDTH
            jrot = cyci_1d([1,nrots], j)
            xsh  = -SHC_INPL_TRSHWDTH
            do while( xsh <= SHC_INPL_TRSHWDTH )
                ysh = -SHC_INPL_TRSHWDTH
                do while( ysh <= SHC_INPL_TRSHWDTH )
                    if( abs(xsh) < SMALL .and. abs(ysh) < SMALL )then
                    else
                        cnt = cnt + 1
                        self%roinds(cnt)   = jrot
                        self%shifts(cnt,1) = xsh
                        self%shifts(cnt,2) = ysh
                    endif
                    ysh = ysh + SHC_INPL_TRSSTEPSZ
                end do
                xsh = xsh + SHC_INPL_TRSSTEPSZ
            end do
        end do
        ! GENERATE THE RANDOM SEARCH ORDER
        call self%rt%reset
        call self%rt%ne_ran_iarr( self%srch_order )
    end subroutine update_grid

    subroutine srch( self, pftcc, ref, iptcl, nrots, prev_rot, rot, shvec )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(shc_inplane),      intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: ref, iptcl, nrots, prev_rot
        integer,                 intent(out)   :: rot
        real,                    intent(out)   :: shvec(2)
        real    :: cc_prev, cc
        integer :: i, ii
        rot     = prev_rot
        shvec   = [0.,0.]
        cc_prev = pftcc%corr(ref, iptcl, prev_rot, shvec)
        call self%update_grid(nrots, prev_rot)
        do i=1,self%maxevals
            ii = self%srch_order(i)
            cc = pftcc%corr(ref, iptcl, self%roinds(ii), self%shifts(ii,:))
            if( cc > cc_prev )then
                rot     = self%roinds(ii)
                shvec   = self%shifts(ii,:)
                cc_prev = cc
                if( i >= self%minevals ) exit
            endif
        end do
    end subroutine srch

    subroutine kill( self )
        class(shc_inplane), intent(out) :: self
        if( self%exists )then
            call self%rt%kill
            deallocate(self%srch_order, self%roinds, self%shifts)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_shc_inplane
