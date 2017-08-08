! Cartesian cross-correlation calculation in Fourier space
module simple_cartft_corrcalc
use simple_defs
use simple_params,    only: params
use simple_build,     only: build
use simple_image,     only: image
use simple_projector, only: projector
implicit none

type :: cartft_corrcalc
    private
    class(params), pointer       :: pp=>null()
    type(projector), allocatable :: refvols(:)
    type(image),     allocatable :: img_refs(:)
    type(image)                  :: img_ctf
    real             :: lp = 0., hp = 0.
    real             :: kv_prev =0., cs_prev =0., fraca_prev =0.
    real             :: dfx_prev=0., dfy_prev=0., angast_prev=0.
    integer          :: nstates   = 0
    character(len=4) :: ctf = 'no'
    logical          :: existence = .false.
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! CTF IMAGE CREATOR
    procedure, private :: create_ctf_image
    ! PROJECTORS
    procedure, private :: project_1
    procedure, private :: project_2
    generic            :: project => project_1, project_2
    ! CORRELATORS
    procedure          :: frc
    procedure, private :: correlate_1
    procedure, private :: correlate_2
    generic            :: correlate => correlate_1, correlate_2
    ! GETTERS
    procedure          :: get_hp
    procedure          :: get_lp
    procedure          :: get_ref
    procedure          :: get_ctf_img
    ! DESTRUCTOR
    procedure          :: kill
end type cartft_corrcalc

contains
    
    !>  \brief  is a constructor
    subroutine new( self, b, p, cline )
        use simple_jiffys,          only: alloc_err
        use simple_hadamard_common, only: preprefvol
        use simple_cmdline,         only: cmdline
        class(cartft_corrcalc), intent(inout) :: self
        class(build),  target,  intent(inout) :: b
        class(params), target,  intent(inout) :: p
        class(cmdline),         intent(inout) :: cline
        integer :: s, alloc_stat
        call self%kill
        ! set constants
        self%nstates = p%nstates
        self%ctf     = p%ctf
        self%lp      = p%lp
        self%hp      = p%hp
        self%pp      => p
        ! allocate reference volumes & one reference image
        allocate(self%refvols(self%nstates), self%img_refs(1), stat=alloc_stat)
        call alloc_err("In: simple_cartft_corrcalc :: new", alloc_stat)
        call self%img_refs(1)%new([self%pp%boxmatch,self%pp%boxmatch,1],self%pp%smpd)
        write(*,'(A)') '>>> PREPARING VOLUMES'
        do s=1,self%nstates
            call preprefvol( b, p, cline, s, doexpand=.false. )
            self%refvols(s) = b%vol
            call self%refvols(s)%expand_cmat
        end do
        self%existence = .true.
    end subroutine new

    !>  \brief  is for creating the CTF image
    subroutine create_ctf_image( self, o )
        use simple_math, only: euclid
        use simple_ori,  only: ori
        use simple_ctf,  only: ctf
        class(cartft_corrcalc), intent(inout) :: self
        class(ori),             intent(inout) :: o
        type(ctf) :: tfun
        real      :: kV, cs, fraca, dfx, dfy, angast, dist
        dfx = o%get('dfx')
        if( self%pp%tfplan%mode .eq. 'astig' )then ! astigmatic CTF
            dfy    = o%get('dfy')
            angast = o%get('angast')
        else if( self%pp%tfplan%mode .eq. 'noastig' )then
            dfy    = dfx
            angast = 0.
        else
            stop 'unsupported ctf mode; create_ctf_image; simple_cartft_corrcalc'
        endif
        kV    = o%get('kv')
        cs    = o%get('cs')
        fraca = o%get('fraca')
        dist = euclid([kV,cs,fraca,dfx,dfy,angast],&
        [self%kv_prev,self%cs_prev,self%fraca_prev,self%dfx_prev,self%dfy_prev,self%angast_prev])
        if( dist < 0.001 )then
            return
        else
            ! CTF parameters have changed, update CTF image
            tfun = ctf(self%pp%smpd, kV, cs, fraca)
            call self%img_ctf%new([self%pp%boxmatch,self%pp%boxmatch,1],self%pp%smpd)
            call tfun%ctf2img(self%img_ctf, dfx, 'ctf', dfy, angast)
            self%kv_prev     = kV 
            self%cs_prev     = cs
            self%fraca_prev  = fraca
            self%dfx_prev    = dfx
            self%dfy_prev    = dfy
            self%angast_prev = angast
        endif
    end subroutine create_ctf_image

    !>  \brief  is for projecting a set
    subroutine project_1( self, os )
        use simple_oris, only: oris
        use simple_ori,  only: ori
        class(cartft_corrcalc), intent(inout) :: self
        class(oris),            intent(inout) :: os
        type(ori) :: o
        integer   :: noris, nrefs, iref
        ! make the container
        noris = os%get_noris()
        nrefs = size(self%img_refs)
        if( noris == nrefs )then
            ! all good
        else
            do iref=1,nrefs
                call self%img_refs(iref)%kill
            end do 
            deallocate(self%img_refs)
            nrefs = noris
            allocate(self%img_refs(nrefs))
            do iref=1,nrefs
                call self%img_refs(iref)%new([self%pp%boxmatch,self%pp%boxmatch,1],self%pp%smpd)
            end do
        endif
        ! project & deal with CTF
        do iref=1,nrefs
            o = os%get_ori(iref)
            call self%project_2(o, iref)            
        end do
    end subroutine project_1

    !>  \brief  is for projecting a single
    subroutine project_2( self, o, iref )
        use simple_ori, only: ori
        class(cartft_corrcalc), intent(inout) :: self
        class(ori),             intent(inout) :: o
        integer,                intent(in)    :: iref
        integer   :: s
        s = nint(o%get('state'))
        call self%refvols(s)%fproject(o, self%img_refs(iref))
        if( self%ctf .ne. 'no' )then
            call self%create_ctf_image(o)
            call self%img_refs(iref)%mul(self%img_ctf)
        endif
    end subroutine project_2

    !>  \brief  for calculating the Fourier Ring Correlation (FRC) betwen the reference
    !!          and the particle
    subroutine frc( self, o, iref, pimg, res, corrs )
        use simple_ori, only: ori
        class(cartft_corrcalc), intent(inout) :: self
        class(ori),             intent(inout) :: o
        integer,                intent(in)    :: iref
        class(image),           intent(inout) :: pimg
        real, allocatable,      intent(out)   :: res(:), corrs(:)
        integer :: s
        s = nint(o%get('state'))
        call self%refvols(s)%fproject(o, self%img_refs(iref))
        if( self%ctf .ne. 'no' )then
            call self%create_ctf_image(o)
            call self%img_refs(iref)%mul(self%img_ctf)
        endif
        call self%img_refs(iref)%fsc(pimg, res, corrs)
    end subroutine frc

    !>  \brief  for calculating all correlations over img_refs (no shifts)
    !!          parameterised over rotational orientations + state
    function correlate_1( self, pimg ) result( cc )
        class(cartft_corrcalc), intent(inout) :: self
        class(image),           intent(inout) :: pimg
        real, allocatable :: cc(:)
        integer :: iref, nrefs
        nrefs = size(self%img_refs)
        allocate(cc(nrefs))
        do iref=1,nrefs
            cc(iref) = self%img_refs(iref)%corr(pimg, self%lp, self%hp)
        end do
    end function correlate_1

    !>  \brief  for calculating the shifted correlation for one reference (iref)
    !!          parameterised over one rotational orientation + state + shift
    function correlate_2( self, pimg, iref, shvec ) result( cc )
        class(cartft_corrcalc), intent(inout) :: self
        class(image),           intent(inout) :: pimg
        integer,                intent(in)    :: iref
        real,                   intent(in)    :: shvec(3)
        real :: cc
        cc = self%img_refs(iref)%corr_shifted(pimg, shvec, self%lp, self%hp)
    end function correlate_2
    
    !>  \brief for getting the current fwd ft CTF img
    function get_ctf_img( self )result( img )
        class(cartft_corrcalc), intent(inout) :: self
        type(image) :: img
        img = self%img_ctf
    end function get_ctf_img

    !>  \brief for getting a reference fwd ft img
    function get_ref( self, ind )result( img )
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: ind
        type(image) :: img
        if( ind<1 .or. ind>size(self%img_refs) )&
            & stop 'Index out of bounds in simple_cartft_corrcalc%get_ref'
        img = self%img_refs( ind )
    end function get_ref

    function get_lp( self )result( lp )
        class(cartft_corrcalc), intent(inout) :: self
        real :: lp
        lp = self%lp
    end function get_lp

    function get_hp( self )result( hp )
        class(cartft_corrcalc), intent(inout) :: self
        real :: hp
        hp = self%hp
    end function get_hp

    !>  \brief  is a destructor
    subroutine kill( self )
        class(cartft_corrcalc), intent(inout) :: self
        integer :: s, iref
        if( self%existence )then
            self%pp => null()
            do s=1,self%nstates
                call self%refvols(s)%kill_expanded
                call self%refvols(s)%kill
            end do
            deallocate(self%refvols)
            if( allocated(self%img_refs) )then
                do iref=1,size(self%img_refs)
                    call self%img_refs(iref)%kill
                end do
                deallocate(self%img_refs)
            endif
            call self%img_ctf%kill
            self%kv_prev     = 0.
            self%cs_prev     = 0.
            self%fraca_prev  = 0.
            self%dfx_prev    = 0.
            self%dfy_prev    = 0.
            self%angast_prev = 0.
            self%lp = 0.
            self%hp = 0.
            self%existence = .false.
        endif
    end subroutine kill

end module simple_cartft_corrcalc
