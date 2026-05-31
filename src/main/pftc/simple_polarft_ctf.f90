!@descr: polarft class submodule for dealing with CTF-related things
submodule (simple_polarft_calc) simple_polarft_ctf
#include "simple_local_flags.inc"
implicit none

contains

    module subroutine create_polar_absctfmats(self, spproj, oritype, pfromto)
        class(polarft_calc),       intent(inout) :: self
        class(sp_project), target, intent(inout) :: spproj
        character(len=*),          intent(in)    :: oritype
        integer, optional,         intent(in)    :: pfromto(2)
        type(ctfparams) :: ctfparms(nthr_glob)
        type(ctf)       :: tfuns(nthr_glob)
        real(sp)        :: spaFreqSq_mat(self%pftsz,self%kfromto(1):self%interpklim)
        real(sp)        :: ang_mat(self%pftsz,self%kfromto(1):self%interpklim), hinv,kinv
        integer         :: i,irot,k,iptcl,ithr,ppfromto(2),ctfmatind,imodel
        logical         :: present_pfromto
        present_pfromto = present(pfromto)
        ppfromto = self%pfromto
        if( present_pfromto ) ppfromto = pfromto
        if( .not. self%ctf_model_audit ) self%nctf_models_seen = 0
        if(.not. self%with_ctf ) return
        if( self%ctf_model_audit )then
            if( .not. allocated(self%ctfparams_ptcls) ) allocate(self%ctfparams_ptcls(1:self%nptcls))
            if( .not. allocated(self%ctfparams_ptcls_set) ) allocate(self%ctfparams_ptcls_set(1:self%nptcls))
            if( .not. allocated(self%ctf_models_seen) ) allocate(self%ctf_models_seen(max(1,self%nptcls)))
            self%ctfparams_ptcls_set = .false.
        endif
        if( self%ctf_scoring_audit )then
            if( .not. allocated(self%ctfparams_scored) ) allocate(self%ctfparams_scored(1:self%nptcls))
            if( .not. allocated(self%ctf_models_scored) ) allocate(self%ctf_models_scored(max(1,self%nptcls)))
            self%ctfparams_scored = .false.
        endif
        !$omp parallel do default(shared) private(irot,k,hinv,kinv) schedule(static) proc_bind(close)
        do irot=1,self%pftsz
            do k=self%kfromto(1),self%interpklim
                hinv = self%polar(1,k,irot) / self%ldim(1)
                kinv = self%polar(2,k,irot) / self%ldim(2)
                spaFreqSq_mat(irot,k) = hinv*hinv+kinv*kinv
                ang_mat(irot,k)       = atan2(self%polar(2,k,irot),self%polar(1,k,irot))
            end do
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(i,iptcl,ctfmatind,ithr,imodel) schedule(static) proc_bind(close)
        do i=ppfromto(1),ppfromto(2)
            if( .not. present_pfromto )then
                iptcl     = i
                ctfmatind = i
            else
                iptcl     = i
                ctfmatind = i - ppfromto(1) + 1
            endif
            if( self%pinds(iptcl) > 0 )then
                ithr           = omp_get_thread_num() + 1
                ctfparms(ithr) = spproj%get_ctfparams(trim(oritype), iptcl)
                tfuns(ithr)    = ctf(ctfparms(ithr)%smpd, ctfparms(ithr)%kv, ctfparms(ithr)%cs, ctfparms(ithr)%fraca)
                call tfuns(ithr)%init(ctfparms(ithr)%dfx, ctfparms(ithr)%dfy, ctfparms(ithr)%angast)
                imodel = self%pinds(ctfmatind)
                if( self%ctf_model_audit )then
                    self%ctfparams_ptcls(imodel)     = ctfparms(ithr)
                    self%ctfparams_ptcls_set(imodel) = .true.
                endif
                if( ctfparms(ithr)%l_phaseplate )then
                    self%ctfmats(:,:,imodel) = abs(tfuns(ithr)%eval(spaFreqSq_mat(:,:), ang_mat(:,:), ctfparms(ithr)%phshift) )
                else
                    self%ctfmats(:,:,imodel) = abs(tfuns(ithr)%eval(spaFreqSq_mat(:,:), ang_mat(:,:)))
                endif
            endif
        end do
        !$omp end parallel do
        if( self%ctf_model_audit ) call count_ctf_models_seen(self%ctfparams_ptcls, self%ctfparams_ptcls_set)
    contains

        subroutine count_ctf_models_seen(ctfmodels_consumed, model_set)
            type(ctfparams), intent(in) :: ctfmodels_consumed(:)
            logical,         intent(in) :: model_set(:)
            integer :: i
            if( size(ctfmodels_consumed) < 1 ) return
            do i = 1,size(ctfmodels_consumed)
                if( .not. model_set(i) ) cycle
                call add_ctf_model_seen(ctfmodels_consumed(i))
            enddo
        end subroutine count_ctf_models_seen

        subroutine add_ctf_model_seen(ctfparms)
            type(ctfparams), intent(in) :: ctfparms
            type(ctfparams), allocatable :: tmp(:)
            integer :: i, oldsz
            do i = 1,self%nctf_models_seen
                if( same_ctf_model(ctfparms, self%ctf_models_seen(i)) ) return
            enddo
            if( self%nctf_models_seen == size(self%ctf_models_seen) )then
                oldsz = size(self%ctf_models_seen)
                allocate(tmp(max(1, 2 * oldsz)))
                tmp(1:oldsz) = self%ctf_models_seen
                call move_alloc(tmp, self%ctf_models_seen)
            endif
            self%nctf_models_seen = self%nctf_models_seen + 1
            self%ctf_models_seen(self%nctf_models_seen) = ctfparms
        end subroutine add_ctf_model_seen

        logical function same_ctf_model(lhs, rhs)
            type(ctfparams), intent(in) :: lhs, rhs
            real, parameter :: CTFTOL = 1.0e-5
            same_ctf_model = (lhs%ctfflag == rhs%ctfflag) .and. &
                (lhs%l_phaseplate .eqv. rhs%l_phaseplate) .and. &
                (abs(lhs%kv    - rhs%kv)    <= CTFTOL) .and. &
                (abs(lhs%cs    - rhs%cs)    <= CTFTOL) .and. &
                (abs(lhs%fraca - rhs%fraca) <= CTFTOL)
        end function same_ctf_model
    end subroutine create_polar_absctfmats

end submodule simple_polarft_ctf
