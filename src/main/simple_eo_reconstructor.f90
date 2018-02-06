! 3D reconstruction of even-odd pairs for FSC estimation
module simple_eo_reconstructor
#include "simple_lib.f08"
use simple_reconstructor, only: reconstructor
use simple_image,         only: image
use simple_params,        only: params
use simple_cmdline,       only: cmdline
use simple_imghead,       only: find_ldim_nptcls
use simple_imgfile,       only: imgfile
use simple_kbinterpol,    only: kbinterpol
use simple_masker,        only: masker
implicit none

public :: eo_reconstructor
private

type :: eo_reconstructor
    private
    type(reconstructor):: even
    type(reconstructor):: odd
    type(reconstructor):: eosum
    type(masker)       :: envmask
    character(len=4)   :: ext
    real               :: res_fsc05          !< target resolution at FSC=0.5
    real               :: res_fsc0143        !< target resolution at FSC=0.143
    real               :: smpd, msk, fny, inner=0., width=10.
    integer            :: box=0, nstates=1, numlen=2
    integer            :: cyc_lims(3,2)  = 0 !< redundant limits
    logical            :: phaseplate = .false.
    logical            :: automsk    = .false.
    logical            :: wiener     = .false.
    logical            :: exists     = .false.
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! SETTERS
    procedure          :: reset_all
    procedure          :: reset_eos
    procedure, private :: reset_eoexp
    procedure, private :: reset_even
    procedure, private :: reset_odd
    procedure          :: reset_sum
    procedure          :: apply_weight
    ! GETTERS
    procedure          :: get_kbwin
    procedure          :: get_res
    ! I/O
    ! writers
    procedure          :: write_eos
    procedure, private :: write_even
    procedure, private :: write_odd
    ! readers
    procedure          :: read_eos
    procedure, private :: read_even
    procedure, private :: read_odd
    ! INTERPOLATION
    procedure, private :: grid_fplane_1
    procedure, private :: grid_fplane_2
    generic            :: grid_fplane => grid_fplane_1, grid_fplane_2
    procedure          :: compress_exp
    procedure          :: expand_exp
    procedure          :: sum_eos !< for merging even and odd into sum
    procedure          :: sum     !< for summing eo_recs obtained by parallel exec
    procedure          :: sampl_dens_correct_eos
    procedure          :: sampl_dens_correct_sum
    ! RECONSTRUCTION
    procedure          :: eorec_distr
    ! DESTRUCTORS
    procedure          :: kill_exp
    procedure          :: kill
end type eo_reconstructor

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine new(self, p )
        class(eo_reconstructor), intent(inout) :: self !< instance
        class(params), target,   intent(in)    :: p    !< parameters object (provides constants)
        logical     :: neg
        call self%kill
        ! set constants
        neg = .false.
        if( p%neg .eq. 'yes' ) neg = .true.
        self%box        = p%box
        self%smpd       = p%smpd
        self%nstates    = p%nstates
        self%inner      = p%inner
        self%width      = p%width
        self%fny        = p%fny
        self%ext        = p%ext
        self%numlen     = p%numlen
        self%msk        = p%msk
        self%automsk    = file_exists(p%mskfile)
        self%phaseplate = p%tfplan%l_phaseplate
        ! create composites
        if( self%automsk )then
            call self%envmask%new([p%box,p%box,p%box], p%smpd)
            call self%envmask%read(p%mskfile)
            call self%envmask%resmask(p)
        endif
        call self%even%new([p%boxpd,p%boxpd,p%boxpd], p%smpd)
        call self%even%alloc_rho(p)
        call self%even%set_ft(.true.)
        call self%odd%new([p%boxpd,p%boxpd,p%boxpd], p%smpd)
        call self%odd%alloc_rho(p)
        call self%odd%set_ft(.true.)
        call self%eosum%new([p%boxpd,p%boxpd,p%boxpd], p%smpd)
        call self%eosum%alloc_rho(p, expand=.false.)
        ! set redundant limits
        self%cyc_lims = self%even%loop_lims(3)
        ! set existence
        self%exists = .true.
    end subroutine new

    ! SETTERS

    !>  \brief  resets all
    subroutine reset_all( self )
        class(eo_reconstructor), intent(inout) :: self
        call self%reset_eos
        call self%reset_eoexp
        call self%reset_sum
    end subroutine reset_all

    !>  \brief  resets the even odd pairs
    subroutine reset_eos( self )
        class(eo_reconstructor), intent(inout) :: self
        call self%even%reset
        call self%odd%reset
    end subroutine reset_eos

    !>  \brief  resets the even odd pairs expanded matrices
    subroutine reset_eoexp( self )
        class(eo_reconstructor), intent(inout) :: self
        call self%even%reset_exp
        call self%odd%reset_exp
    end subroutine reset_eoexp

    !>  \brief  resets the even
    subroutine reset_even( self )
        class(eo_reconstructor), intent(inout) :: self
        call self%even%reset
    end subroutine reset_even

    !>  \brief  resets the odd
    subroutine reset_odd( self )
        class(eo_reconstructor), intent(inout) :: self
        call self%odd%reset
    end subroutine reset_odd

    !>  \brief  resets the sum
    subroutine reset_sum( self )
        class(eo_reconstructor), intent(inout) :: self
        call self%eosum%reset
    end subroutine reset_sum

    subroutine apply_weight( self, w )
        class(eo_reconstructor), intent(inout) :: self
        real,                    intent(in)    :: w
        call self%even%apply_weight(w)
        call self%odd%apply_weight(w)
    end subroutine apply_weight

    ! GETTERS

    !>  \brief  return the window functions used by eo_reconstructor
    function get_kbwin( self ) result( wf )
        use simple_kbinterpol, only: kbinterpol
        class(eo_reconstructor), intent(inout) :: self
        type(kbinterpol) :: wf
        wf = self%even%get_kbwin()
    end function get_kbwin

    !> \brief  for getting the resolution
    !> \param res_fsc05  target resolution a FSC=0.5
    !> \param res_fsc0143  target resolution a FSC=0.143
    subroutine get_res( self, res_fsc05, res_fsc0143 )
        class(eo_reconstructor), intent(in)  :: self !< instance
        real,                    intent(out) :: res_fsc05, res_fsc0143
        res_fsc0143 = self%res_fsc0143
        res_fsc05   = self%res_fsc05
    end subroutine get_res

    ! I/O

    !>  \brief  write the even and odd reconstructions
    subroutine write_eos( self, fbody )
        class(eo_reconstructor), intent(inout) :: self
        character(len=*),        intent(in)    :: fbody !< filename
        call self%write_even(fbody)
        call self%write_odd(fbody)
    end subroutine write_eos

    !>  \brief  write the even reconstruction
    subroutine write_even( self, fbody )
        class(eo_reconstructor), intent(inout) :: self
        character(len=*),        intent(in)    :: fbody
        call self%even%write(trim(adjustl(fbody))//'_even'//self%ext, del_if_exists=.true.)
        call self%even%write_rho(trim('rho_'//trim(adjustl(fbody))//'_even'//self%ext))
    end subroutine write_even

    !>  \brief  write the odd reconstruction
    subroutine write_odd( self, fbody )
        class(eo_reconstructor), intent(inout) :: self
        character(len=*),        intent(in)    :: fbody
        call self%odd%write(trim(adjustl(fbody))//'_odd'//self%ext, del_if_exists=.true.)
        call self%odd%write_rho('rho_'//trim(adjustl(fbody))//'_odd'//self%ext)
    end subroutine write_odd

    !>  \brief read the even and odd reconstructions
    subroutine read_eos( self, fbody )
        class(eo_reconstructor), intent(inout) :: self
        character(len=*),        intent(in)    :: fbody
        call self%read_even(fbody)
        call self%read_odd(fbody)
    end subroutine read_eos

    !>  \brief  read the even reconstruction
    subroutine read_even( self, fbody )
        class(eo_reconstructor), intent(inout) :: self
        character(len=*),        intent(in)    :: fbody
        character(len=STDLEN)                  :: even_vol, even_rho
        logical                                :: here(2)
        even_vol = trim(adjustl(fbody))//'_even'//self%ext
        even_rho = 'rho_'//trim(adjustl(fbody))//'_even'//self%ext
        here(1)= file_exists(even_vol)
        here(2)= file_exists(even_rho)
        if( all(here) )then
            call self%even%read(even_vol)
            call self%even%read_rho(even_rho)
        else
            call self%reset_even
        endif
    end subroutine read_even

    !>  \brief  read the odd reconstruction
    subroutine read_odd( self, fbody )
        class(eo_reconstructor), intent(inout) :: self
        character(len=*),        intent(in)    :: fbody
        character(len=STDLEN)                  :: odd_vol, odd_rho
        logical                                :: here(2)
        odd_vol = trim(adjustl(fbody))//'_odd'//self%ext
        odd_rho = 'rho_'//trim(adjustl(fbody))//'_odd'//self%ext
        here(1)= file_exists(odd_vol)
        here(2)= file_exists(odd_rho)
        if( all(here) )then
            call self%odd%read(odd_vol)
            call self%odd%read_rho(odd_rho)
        else
            call self%reset_odd
        endif
    end subroutine read_odd

    ! INTERPOLATION

    !> \brief  for gridding a Fourier plane
    subroutine grid_fplane_1( self, se, o, fpl, eo, pwght )
        use simple_ori, only: ori
        use simple_sym, only: sym
        class(eo_reconstructor), intent(inout) :: self  !< instance
        class(sym),              intent(inout) :: se    !< symmetry elements
        class(ori),              intent(inout) :: o     !< orientation
        class(image),            intent(inout) :: fpl   !< Fourier plane
        integer,                 intent(in)    :: eo    !< eo flag
        real,                    intent(in)    :: pwght !< external particle weight (affects both fplane and rho)
        select case(eo)
            case(-1,0)
                call self%even%insert_fplane(se, o, fpl, pwght)
            case(1)
                call self%odd%insert_fplane(se, o, fpl, pwght)
            case DEFAULT
                stop 'unsupported eo flag; eo_reconstructor :: grid_fplane'
        end select
    end subroutine grid_fplane_1

    subroutine grid_fplane_2( self, se, os, fpl, eo, pwght, state )
        use simple_oris, only: oris
        use simple_sym,  only: sym
        class(eo_reconstructor), intent(inout) :: self  !< instance
        class(sym),              intent(inout) :: se    !< symmetry elements
        class(oris),             intent(inout) :: os    !< orientation
        class(image),            intent(inout) :: fpl   !< Fourier plane
        integer,                 intent(in)    :: eo    !< eo flag
        real,                    intent(in)    :: pwght !< external particle weight (affects both fplane and rho)
        integer, optional,       intent(in)    :: state !< state flag
        select case(eo)
            case(-1,0)
                call self%even%insert_fplane(se, os, fpl, pwght, state)
            case(1)
                call self%odd%insert_fplane(se, os, fpl, pwght, state)
            case DEFAULT
                stop 'unsupported eo flag; eo_reconstructor :: grid_fplane'
        end select
    end subroutine grid_fplane_2

    !> \brief  for summing the even odd pairs, resulting sum in self%even
    subroutine sum_eos( self )
        class(eo_reconstructor), intent(inout) :: self !< instance
        call self%eosum%reset
        call self%eosum%sum(self%even)
        call self%eosum%sum(self%odd)
    end subroutine sum_eos

    !> \brief  for summing reconstructors generated by parallel execution
    subroutine sum( self, self_in )
         class(eo_reconstructor), intent(inout) :: self
         class(eo_reconstructor), intent(in)    :: self_in
         call self%even%sum(self_in%even)
         call self%odd%sum(self_in%odd)
    end subroutine sum

    !>  \brief compress e/o
    subroutine compress_exp( self )
        class(eo_reconstructor), intent(inout) :: self
        call self%even%compress_exp
        call self%odd%compress_exp
    end subroutine compress_exp

    !>  \brief expand e/o
    subroutine expand_exp( self )
        class(eo_reconstructor), intent(inout) :: self
        call self%even%expand_exp
        call self%odd%expand_exp
    end subroutine expand_exp

    !> \brief  for sampling density correction of the eo pairs
    subroutine sampl_dens_correct_eos( self, state, fname_even, fname_odd, resmskname, find4eoavg )
        use simple_masker,  only: masker
        class(eo_reconstructor), intent(inout) :: self                  !< instance
        integer,                 intent(in)    :: state                 !< state
        character(len=*),        intent(in)    :: fname_even, fname_odd !< even/odd filenames
        character(len=*),        intent(in)    :: resmskname            !< resolution mask name
        integer,                 intent(out)   :: find4eoavg            !< Fourier index for eo averaging
        real, allocatable :: res(:), corrs(:)
        type(image)       :: even, odd
        integer           :: j, find_plate
        ! make clipped volumes
        call even%new([self%box,self%box,self%box],self%smpd)
        call odd%new([self%box,self%box,self%box],self%smpd)
        ! correct for the uneven sampling density
        call self%even%sampl_dens_correct(maxits=0)
        call self%odd%sampl_dens_correct(maxits=0)
        ! reverse FT
        call self%even%bwd_ft
        call self%odd%bwd_ft
        ! clip
        call self%even%clip(even)
        call self%odd%clip(odd)
        ! write unnormalised unmasked even/odd volumes
        call even%write(trim(fname_even), del_if_exists=.true.)
        call odd%write(trim(fname_odd),   del_if_exists=.true.)
        ! always normalise before masking
        call even%norm
        call odd%norm
        if( self%automsk )then
            call even%zero_background
            call odd%zero_background
            call even%mul(self%envmask)
            call odd%mul(self%envmask)
            call self%envmask%write(resmskname)
        else
            ! spherical masking
            if( self%inner > 1. )then
                call even%mask(self%msk, 'soft', inner=self%inner, width=self%width)
                call odd%mask(self%msk, 'soft', inner=self%inner, width=self%width)
            else
                call even%mask(self%msk, 'soft')
                call odd%mask(self%msk, 'soft')
            endif
        endif
        ! forward FT
        call even%fwd_ft
        call odd%fwd_ft
        ! calculate FSC
        res = even%get_res()
        allocate(corrs(even%get_filtsz()))
        call even%fsc(odd, corrs)
        find_plate = 0
        if( self%phaseplate ) call phaseplate_correct_fsc(corrs, find_plate)
        do j=1,size(res)
           write(*,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
        end do
        ! save, get & print resolution
        call arr2file(corrs, 'fsc_state'//int2str_pad(state,2)//'.bin')
        call get_resolution(corrs, res, self%res_fsc05, self%res_fsc0143)
        self%res_fsc05   = max(self%res_fsc05,self%fny)
        self%res_fsc0143 = max(self%res_fsc0143,self%fny)
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', self%res_fsc05
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', self%res_fsc0143
        ! Fourier index for eo averaging
        find4eoavg = max(K4EOAVGLB,get_lplim_at_corr(corrs, FSC4EOAVG3D))
        find4eoavg = max(find4eoavg, find_plate)
        deallocate(corrs, res)
        call even%kill
        call odd%kill
    end subroutine sampl_dens_correct_eos

    !> \brief  for sampling density correction, antialiasing, bwd_ft & normalization of the sum
    subroutine sampl_dens_correct_sum( self, reference )
        class(eo_reconstructor), intent(inout) :: self      !< instance
        class(image),            intent(inout) :: reference !< reference volume
        write(*,'(A)') '>>> SAMPLING DENSITY (RHO) CORRECTION & WIENER NORMALIZATION'
        call reference%set_ft(.false.)
        call self%eosum%sampl_dens_correct
        call self%eosum%bwd_ft
        call self%eosum%clip(reference)
    end subroutine sampl_dens_correct_sum

    ! RECONSTRUCTION

    !> \brief  for distributed reconstruction of even/odd maps
    subroutine eorec_distr( self, p, o, se, state, fbody )
        use simple_oris,       only: oris
        use simple_sym,        only: sym
        use simple_params,     only: params
        use simple_prep4cgrid, only: prep4cgrid
        class(eo_reconstructor),    intent(inout) :: self   !< object
        class(params),              intent(in)    :: p      !< parameters
        class(oris),                intent(inout) :: o      !< orientations
        class(sym),                 intent(inout) :: se     !< symmetry element
        integer,                    intent(in)    :: state  !< state to reconstruct
        character(len=*), optional, intent(in)    :: fbody  !< body of output file
        type(kbinterpol)  :: wf
        type(image)       :: img, img_pad
        type(prep4cgrid)  :: gridprep
        real              :: skewness
        integer           :: statecnt(p%nstates), i, cnt, state_here, state_glob
        ! stash global state index
        state_glob = state
        ! make the images
        call img%new([p%box,p%box,1],p%smpd)
        call img_pad%new([p%boxpd,p%boxpd,1],p%smpd)
        ! make the gridding prepper
        wf = self%even%get_kbwin()
        call gridprep%new(img, wf, [p%boxpd,p%boxpd,1])
        ! population balancing logics
        if( p%balance > 0 )then
            call o%balance( p%balance, NSPACE_BALANCE, p%nsym, p%eullims, skewness )
            write(*,'(A,F8.2)') '>>> PROJECTION DISTRIBUTION SKEWNESS(%):', 100. * skewness
        else
            call o%set_all2single('state_balance', 1.0)
        endif
        ! zero the Fourier volumes and rhos
        call self%reset_all
        call self%reset_eoexp
        write(*,'(A)') '>>> KAISER-BESSEL INTERPOLATION'
        statecnt = 0
        cnt      = 0
        do i=1,p%nptcls
            call progress(i, p%nptcls)
            if( i <= p%top .and. i >= p%fromp )then
                cnt = cnt + 1
                state_here = nint(o%get(i,'state'))
                if( state_here > 0 .and. (state_here == state ) )then
                    statecnt(state) = statecnt(state) + 1
                    call rec_dens
                endif
            endif
        end do
        ! undo fourier components expansion
        call self%compress_exp
        ! density correction & output
        if( p%l_distr_exec )then
            if( present(fbody) )then
                call self%write_eos(fbody//int2str_pad(state,2)//'_part'//int2str_pad(p%part,self%numlen))
            else
                call self%write_eos('recvol_state'//int2str_pad(state,2)//'_part'//int2str_pad(p%part,self%numlen))
            endif
        endif
        call img%kill
        call img_pad%kill
        ! report how many particles were used to reconstruct each state
        if( p%nstates > 1 )then
            write(*,'(a,1x,i3,1x,a,1x,i6)') '>>> NR OF PARTICLES INCLUDED IN STATE:', state, 'WAS:', statecnt(state)
        endif

        contains

            !> \brief  the density reconstruction functionality
            subroutine rec_dens
                use simple_ori, only: ori
                character(len=:), allocatable :: stkname
                type(ori) :: orientation
                integer   :: state, state_balance, ind, eo
                real      :: pw
                state         = nint(o%get(i, 'state'))
                state_balance = nint(o%get(i, 'state_balance'))
                if( state == 0 .or. state_balance == 0 ) return
                pw = 1.
                if( p%frac < 0.99 ) pw = o%get(i, 'w')
                if( pw > 0. )then
                    orientation = o%get_ori(i)
                    eo          = nint(orientation%get('eo'))
                    if( p%l_stktab_input )then
                        call p%stkhandle%get_stkname_and_ind(i, stkname, ind)
                    else
                        if( p%l_distr_exec )then
                            ind = cnt
                            allocate(stkname, source=trim(p%stk_part))
                        else
                            ind = i
                            allocate(stkname, source=trim(p%stk))
                        endif
                    endif
                    call img%read(stkname, ind)
                    call gridprep%prep(img, img_pad)
                    call self%grid_fplane(se, orientation, img_pad, eo, pw)
                 endif
            end subroutine rec_dens

    end subroutine eorec_distr

    ! DESTRUCTORS

    !>  \brief  is the expanded destructor
    subroutine kill_exp( self )
        class(eo_reconstructor), intent(inout) :: self !< instance
        if( self%exists )then
            call self%even%dealloc_exp
            call self%odd%dealloc_exp
            call self%eosum%dealloc_exp
        endif
    end subroutine kill_exp

    !>  \brief  is a destructor
    subroutine kill( self )
        class(eo_reconstructor), intent(inout) :: self !< instance
        if( self%exists )then
            ! kill composites
            call self%envmask%kill
            call self%even%dealloc_rho
            call self%even%kill
            call self%odd%dealloc_rho
            call self%odd%kill
            call self%eosum%dealloc_rho
            call self%eosum%kill
            ! set existence
            self%exists = .false.
        endif
    end subroutine kill

end module simple_eo_reconstructor
