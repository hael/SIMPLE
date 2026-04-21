!@descr: 2D strategy for probabilistic class assignment (precomputed by prob_align2D/prob_tab2D)
module simple_strategy2D_prob
use simple_pftc_srch_api
use simple_eul_prob_tab,     only: eulprob_corr_switch
use simple_parameters,       only: parameters
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_spec
use simple_oris,             only: oris
use simple_builder,          only: builder
implicit none

public :: strategy2D_prob
private

#include "simple_local_flags.inc"

type, extends(strategy2D) :: strategy2D_prob
contains
    procedure :: new  => new_prob
    procedure :: srch => srch_prob
    procedure :: kill => kill_prob
end type strategy2D_prob

contains

    subroutine new_prob( self, params, spec, build )
        class(strategy2D_prob), intent(inout) :: self
        class(parameters),      intent(in)    :: params
        class(strategy2D_spec), intent(inout) :: spec
        class(builder),         intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_prob

    subroutine srch_prob( self, os )
        class(strategy2D_prob), intent(inout) :: self
        class(oris),            intent(inout) :: os
        integer :: iptcl_map, icls, inpl
        real    :: corr
        if( os%get_state(self%s%iptcl) > 0 )then
            if( .not. associated(self%spec%eulprob_obj_part2D) ) THROW_HARD('strategy2D_prob requires eulprob_obj_part2D')
            call self%s%prep4srch(os)
            self%s%nrefs_eval = self%s%nrefs
            iptcl_map = self%s%iptcl_map
            icls      = self%spec%eulprob_obj_part2D%assgn_map(iptcl_map)%icls
            inpl      = self%spec%eulprob_obj_part2D%assgn_map(iptcl_map)%inpl
            corr      = eulprob_corr_switch(self%spec%eulprob_obj_part2D%assgn_map(iptcl_map)%dist, self%s%p_ptr%cc_objfun)
            self%s%best_class = icls
            self%s%best_rot   = inpl
            self%s%best_corr  = corr
            self%s%best_shvec = [0.,0.]
            if( self%s%p_ptr%l_doshift .and. self%spec%eulprob_obj_part2D%assgn_map(iptcl_map)%has_sh )then
                self%s%best_shvec = [self%spec%eulprob_obj_part2D%assgn_map(iptcl_map)%x, &
                                  self%spec%eulprob_obj_part2D%assgn_map(iptcl_map)%y]
            endif
            call self%s%store_solution(self%s%best_class, self%s%best_rot, self%s%best_corr)
            call self%s%assign_ori(os)
        else
            call os%reject(self%s%iptcl)
        endif
    end subroutine srch_prob

    subroutine kill_prob( self )
        class(strategy2D_prob), intent(inout) :: self
        call self%s%kill
    end subroutine kill_prob

end module simple_strategy2D_prob
