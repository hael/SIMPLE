!@descr: the abstract data type for sigma2 used when objfun=euclid
module simple_euclid_sigma2
use, intrinsic :: iso_fortran_env, only: int64, real32
use simple_core_module_api
use simple_polarft_calc,   only: polarft_calc
use simple_cartft_calc,    only: cartft_calc
use simple_fplane_grid,    only: fplane_grid
use simple_parameters,     only: parameters
use simple_syslib,         only: simple_atomic_replace
use simple_sigma2_state,   only: sigma2_state_range_path, sigma2_state_next_generation
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_read_header, &
    &sigma2_state_read_global, sigma2_state_read_particles_into, sigma2_state_write_local_range, &
    &sigma2_state_read_local_range, SIGMA2_STATE_COMMITTED
implicit none
private

#include "simple_local_flags.inc"

public :: euclid_sigma2, sigma2_commit_iter, sigma2_half

type euclid_sigma2
    private
    class(parameters),    pointer :: p_ptr => null()
    real,    allocatable, public  :: sigma2_noise(:,:)      !< the two global spectra (shell, half): column sigma2_half(eo)
    real,    allocatable          :: sigma2_part(:,:)       !< the residual records of this part (matchers only)
    integer                       :: fromp
    integer                       :: top
    integer                       :: kfromto(2) = 0
    type(string)                  :: binfname
    logical                       :: exists     = .false.
    logical                       :: replaces_range = .false. !< started from the open transaction's range (read_pending_range)
    ! the reconstruction-plane grid of the spectra (shell table, resampled inverse spectra), memoised by
    ! prep_fplane_grid for one padded box and lowest shell; dropped whenever the spectra are (re)loaded
    type(fplane_grid), public     :: fplane
contains
    ! constructor
    procedure, private :: new_pftc
    procedure, private :: new_cart
    generic            :: new => new_pftc, new_cart
    ! utils
    procedure          :: get_kfromto
    procedure          :: set_kfromto
    ! I/O
    procedure          :: read_part
    procedure          :: read_pending_range
    procedure          :: read_global
    procedure          :: prep_fplane_grid
    procedure, private :: calc_sigma2_pftc
    procedure, private :: calc_sigma2_cart
    procedure, private :: calc_sigma2_cart_residual
    generic            :: calc_sigma2 => calc_sigma2_pftc, calc_sigma2_cart, calc_sigma2_cart_residual
    procedure, private :: check_record
    procedure          :: get_sigma2_part
    procedure          :: write_sigma2
    ! destructor
    procedure          :: kill
end type euclid_sigma2

contains

    pure integer function sigma2_commit_iter( matcher_iter, matcher_completed ) result( commit_iter )
        integer, intent(in) :: matcher_iter
        logical, intent(in) :: matcher_completed
        commit_iter = matcher_iter
        if( matcher_completed ) commit_iter = commit_iter + 1
    end function sigma2_commit_iter

    !> The column of a particle's half (eo 0 even, 1 odd) in the two global spectra; any other
    !! assignment is refused, as the reduction refuses it
    integer function sigma2_half( eo ) result( half )
        integer, intent(in) :: eo
        if( eo < 0 .or. eo > 1 ) THROW_HARD('particle without an even/odd assignment (eo='//int2str(eo)//') has no sigma2 spectrum')
        half = eo + 1
    end function sigma2_half

    !> Polar pass: the noise table is shared with the polar calculator.
    subroutine new_pftc( self, params, pftc, binfname, box )
        class(euclid_sigma2), target, intent(inout) :: self
        class(parameters),    target, intent(in)    :: params
        class(polarft_calc),          intent(inout) :: pftc
        class(string),                intent(in)    :: binfname
        integer,                      intent(in)    :: box
        call self%new_cart(params, binfname, box)
        call pftc%assign_sigma2_noise(self%sigma2_noise)
    end subroutine new_pftc

    !> The two global spectra of a pass: a Cartesian pass, which constructs no polar calculator,
    !! uses them alone; the polar constructor shares them with the calculator.
    subroutine new_cart( self, params, binfname, box )
        class(euclid_sigma2), target, intent(inout) :: self
        class(parameters),    target, intent(in)    :: params
        class(string),                intent(in)    :: binfname
        integer,                      intent(in)    :: box

        call self%kill
        self%p_ptr => params
        self%kfromto = [1, fdim(box)-1]
        allocate(self%sigma2_noise(self%kfromto(1):self%kfromto(2),2), source=0.)
        self%binfname = binfname
        self%fromp = self%p_ptr%fromp
        self%top = self%p_ptr%top
        self%exists = .true.
    end subroutine new_cart

    pure function get_kfromto( self )result( kfromto )
        class(euclid_sigma2), intent(in) :: self
        integer :: kfromto(2)
        kfromto = self%kfromto
    end function get_kfromto

    !>  Set the reconstruction/alignment band. Used by workflows (e.g. flex_pca)
    !!  that reconstruct at a cropped box without loaded sigma2 spectra, so the
    !!  reconstruction backend uses the correct box_crop spectral range.
    subroutine set_kfromto( self, kfromto )
        class(euclid_sigma2), intent(inout) :: self
        integer,              intent(in)    :: kfromto(2)
        self%kfromto = kfromto
    end subroutine set_kfromto

    ! I/O

    !> The residual records of this part's rows, which a matcher updates and writes back: read straight
    !! into the record array (shells kfromto, rows fromp:top); the reader refuses another shell count
    subroutine read_part( self )
        class(euclid_sigma2), intent(inout) :: self
        integer :: status
        character(len=STDLEN) :: message
        if( allocated(self%sigma2_part) ) deallocate(self%sigma2_part)
        allocate(self%sigma2_part(self%kfromto(1):self%kfromto(2),self%fromp:self%top))
        call sigma2_state_read_particles_into(self%binfname%to_char(), self%fromp, self%sigma2_part, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
    end subroutine read_part

    !> Replace this part's rows by the range an earlier pass of the open transaction wrote for it
    !! (the discrete pass a polish follows), so a particle the later pass leaves untouched keeps
    !! the residual of its current pose rather than the committed generation's.
    subroutine read_pending_range( self )
        class(euclid_sigma2), intent(inout) :: self
        type(string) :: range_path
        real(real32), allocatable :: spectra(:,:)
        integer(int64) :: next_gen, generation, layout_digest
        integer :: first_row, last_row, kfrom, kto, status
        character(len=STDLEN) :: message
        if( .not. allocated(self%sigma2_part) ) THROW_HARD('particle sigma2 storage is not allocated')
        call sigma2_state_next_generation(self%binfname%to_char(), next_gen, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        range_path = sigma2_state_range_path(self%binfname%to_char(), next_gen, self%p_ptr%part, self%p_ptr%numlen)
        call sigma2_state_read_local_range(range_path%to_char(), generation, layout_digest, first_row, last_row, &
            &kfrom, kto, spectra, status, message)
        if( status /= 0 ) THROW_HARD('no sigma2 range of the open transaction for this part: '//trim(message))
        if( generation /= next_gen ) THROW_HARD('pending sigma2 range belongs to another transaction')
        if( first_row /= self%fromp .or. last_row /= self%top ) THROW_HARD('pending sigma2 range does not cover this part')
        if( kfrom /= self%kfromto(1) .or. kto /= self%kfromto(2) ) THROW_HARD('pending sigma2 range has incompatible shell bounds')
        self%sigma2_part(:,:) = real(spectra)
        self%replaces_range   = .true.
        deallocate(spectra)
        call range_path%kill
    end subroutine read_pending_range

    !> The two global spectra of the committed state, one per half
    subroutine read_global( self )
        class(euclid_sigma2), intent(inout) :: self
        real(real32), allocatable :: spectra(:,:)
        integer               :: status
        character(len=STDLEN) :: message
        if( .not.associated(self%p_ptr) )then
            THROW_HARD('euclid_sigma2: params pointer is not set')
        endif
        call sigma2_state_read_global(self%binfname%to_char(), spectra, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        if( size(spectra,1) /= self%kfromto(2)-self%kfromto(1)+1 ) &
            &THROW_HARD('canonical global sigma2 has incompatible shell bounds')
        self%sigma2_noise(:,:) = real(spectra)
        deallocate(spectra)
        call self%fplane%kill
    end subroutine read_global

    !> The reconstruction-plane grid of a padded box box_croppd with the two spectra attached from shell
    !! kfrom (fplane): built once and kept while the box, the shell and the spectra stay the same
    subroutine prep_fplane_grid( self, box_croppd, kfrom )
        class(euclid_sigma2), intent(inout) :: self
        integer,              intent(in)    :: box_croppd, kfrom
        if( .not. allocated(self%sigma2_noise) ) THROW_HARD('euclid_sigma2 holds no spectra; prep_fplane_grid')
        if( self%fplane%matches(box_croppd, kfrom) ) return
        call self%fplane%new(box_croppd)
        call self%fplane%set_sigma2(lbound(self%sigma2_noise,1), self%sigma2_noise, kfrom)
    end subroutine prep_fplane_grid

    !>  Calculates and updates sigma2 within search resolution range
    !> The sigma2 contribution of particle iptcl at its assigned polar orientation o.
    subroutine calc_sigma2_pftc( self, pftc, iptcl, o, refkind )
        class(euclid_sigma2), intent(inout) :: self
        class(polarft_calc),  intent(inout) :: pftc
        integer,              intent(in)    :: iptcl
        class(ori),           intent(in)    :: o
        character(len=*),     intent(in)    :: refkind ! 'proj' or 'class'
        integer :: iref, irot, kfromto(2)
        real    :: shvec(2)
        if( .not.associated(self%p_ptr) )then
            THROW_HARD('euclid_sigma2: params pointer is not set')
        endif
        if ( o%isstatezero() ) return
        call self%check_record(iptcl)
        kfromto = self%p_ptr%kfromto
        shvec = o%get_2Dshift()
        iref  = nint(o%get(trim(refkind)))
        if( trim(refkind) == 'proj' .and. self%p_ptr%nstates > 1 )then
            iref = (o%get_state() - 1) * self%p_ptr%nspace + iref
        endif
        irot  = pftc%get_roind(360. - o%e3get())
        ! the residual goes straight into the particle's record over the pass's band
        call pftc%gen_sigma_contrib(iref, iptcl, shvec, irot, self%sigma2_part(kfromto(1):kfromto(2),iptcl))
    end subroutine calc_sigma2_pftc

    !> The sigma2 contribution of particle iptcl at the committed Cartesian pose, from particle
    !! slot islot of the Cartesian calculator (prepared for objfun=euclid; sigma2 belongs to the
    !! Euclidean objective only, C5): the mean squared residual over two per sample per shell
    !! of the slot's shell range (cartft_calc%sigma_contribution). rotmat and shift (cropped
    !! pixels, the model shift relative to the slot's observation) are the committed pose;
    !! state and iseven select the reference. THROW_HARD when the slot's shell range is not the
    !! pass's band.
    subroutine calc_sigma2_cart( self, cftc, islot, iptcl, state, iseven, rotmat, shift )
        class(euclid_sigma2), intent(inout) :: self
        class(cartft_calc),   intent(in)    :: cftc
        integer,              intent(in)    :: islot, iptcl, state
        logical,              intent(in)    :: iseven
        real(dp),             intent(in)    :: rotmat(3,3), shift(2)
        real, allocatable :: sigma_contrib(:)
        if( .not.associated(self%p_ptr) )then
            THROW_HARD('euclid_sigma2: params pointer is not set')
        endif
        if( .not. cftc%ptcl_is_valid(islot) ) return
        ! residual only: the diagnostic outputs of sigma_contribution are test instruments
        call cftc%sigma_contribution(state, iseven, islot, rotmat, shift, sigma_contrib)
        if( lbound(sigma_contrib,1) /= self%p_ptr%kfromto(1) .or. ubound(sigma_contrib,1) /= self%p_ptr%kfromto(2) ) &
            &THROW_HARD('Cartesian sigma contribution does not span the band of the pass')
        call self%check_record(iptcl)
        self%sigma2_part(self%p_ptr%kfromto(1):self%p_ptr%kfromto(2),iptcl) = sigma_contrib
    end subroutine calc_sigma2_cart

    !> calc_sigma2_cart from the shell residual of the evaluation that chose the committed pose
    !! (cartft_pose_opt%get_shell_residual), without a pass of its own
    subroutine calc_sigma2_cart_residual( self, cftc, islot, iptcl, shell_residual )
        class(euclid_sigma2), intent(inout) :: self
        class(cartft_calc),   intent(in)    :: cftc
        integer,              intent(in)    :: islot, iptcl
        real(dp),             intent(in)    :: shell_residual(:)
        real, allocatable :: sigma_contrib(:)
        if( .not.associated(self%p_ptr) )then
            THROW_HARD('euclid_sigma2: params pointer is not set')
        endif
        if( .not. cftc%ptcl_is_valid(islot) ) return
        call cftc%sigma_from_shell_residual(islot, shell_residual, sigma_contrib)
        if( lbound(sigma_contrib,1) /= self%p_ptr%kfromto(1) .or. ubound(sigma_contrib,1) /= self%p_ptr%kfromto(2) ) &
            &THROW_HARD('Cartesian sigma contribution does not span the band of the pass')
        call self%check_record(iptcl)
        self%sigma2_part(self%p_ptr%kfromto(1):self%p_ptr%kfromto(2),iptcl) = sigma_contrib
    end subroutine calc_sigma2_cart_residual

    !> The particle's sigma2 contribution of this part over the pass's band (as calc_sigma2 stored it).
    function get_sigma2_part( self, iptcl ) result( sigma2 )
        class(euclid_sigma2), intent(in) :: self
        integer,              intent(in) :: iptcl
        real, allocatable :: sigma2(:)
        if( .not. allocated(self%sigma2_part) ) THROW_HARD('particle sigma2 storage is not allocated')
        if( iptcl < lbound(self%sigma2_part,2) .or. iptcl > ubound(self%sigma2_part,2) ) &
            &THROW_HARD('particle index is outside sigma2 storage')
        sigma2 = self%sigma2_part(self%p_ptr%kfromto(1):self%p_ptr%kfromto(2),iptcl)
    end function get_sigma2_part

    !> A residual can be stored for particle iptcl: the records are read and the row is this part's.
    subroutine check_record( self, iptcl )
        class(euclid_sigma2), intent(in) :: self
        integer,              intent(in) :: iptcl
        if( .not. allocated(self%sigma2_part) ) THROW_HARD('particle sigma2 storage is not allocated')
        if( iptcl < lbound(self%sigma2_part,2) .or. iptcl > ubound(self%sigma2_part,2) ) &
            &THROW_HARD('particle index is outside sigma2 storage')
    end subroutine check_record

    subroutine write_sigma2( self )
        class(euclid_sigma2), intent(inout) :: self
        type(sigma2_state_header) :: header
        type(string) :: range_path, staged_path
        integer :: status
        character(len=STDLEN) :: message
        ! the range carries the identity of the committed store and the generation the master's
        ! commit will make, so nothing of another transaction can be merged into it
        call sigma2_state_read_header(self%binfname%to_char(), header, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        if( header%state /= SIGMA2_STATE_COMMITTED ) THROW_HARD('canonical sigma2 update source is not committed')
        header%generation = header%generation + 1_int64
        range_path = sigma2_state_range_path(self%binfname%to_char(), header%generation, self%p_ptr%part, self%p_ptr%numlen)
        if( self%replaces_range )then
            ! the range this pass started from stays until its successor is complete
            staged_path = range_path%to_char()//'.tmp'
            call del_file(staged_path)
            call sigma2_state_write_local_range(staged_path%to_char(), header%generation, header%layout_digest, &
                &self%fromp, self%sigma2_part, self%kfromto(1), self%kfromto(2), status, message)
            if( status /= 0 ) THROW_HARD(trim(message))
            call simple_atomic_replace(staged_path, range_path, status)
            if( status /= 0 ) THROW_HARD('cannot replace the pending sigma2 range')
            call staged_path%kill
        else
            call sigma2_state_write_local_range(range_path%to_char(), header%generation, header%layout_digest, &
                &self%fromp, self%sigma2_part, self%kfromto(1), self%kfromto(2), status, message)
            if( status /= 0 ) THROW_HARD(trim(message))
        endif
        call range_path%kill
    end subroutine write_sigma2

    ! Destructor

    subroutine kill( self )
        class(euclid_sigma2), intent(inout) :: self
        if( self%exists )then
            if(allocated(self%sigma2_noise))  deallocate(self%sigma2_noise)
            if( allocated(self%sigma2_part) ) deallocate(self%sigma2_part)
            self%kfromto     = 0
            self%fromp       = -1
            self%top         = -1
            self%exists      = .false.
        endif
        call self%fplane%kill
        self%replaces_range = .false.
        self%p_ptr => null()
    end subroutine kill

end module simple_euclid_sigma2
