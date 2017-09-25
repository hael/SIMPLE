! 3D reconstruction of even-odd pairs for FSC estimation
module simple_eo_reconstructor
use simple_defs           ! use all in there
use simple_reconstructor, only: reconstructor
use simple_image,         only: image
use simple_params,        only: params
use simple_cmdline,       only: cmdline
implicit none

public :: eo_reconstructor
private

type :: eo_reconstructor
    private
    type(reconstructor) :: even
    type(reconstructor) :: odd
    type(reconstructor) :: eosum
    type(image)         :: envmask
    character(len=4)    :: ext
    real                :: fsc05      !< target resolution at FSC=0.5
    real                :: fsc0143    !< target resolution at FSC=0.143
    real                :: smpd, msk, fny, inner=0., width=10.
    integer             :: box=0, nstates=1, numlen=2, lfny=0
    logical             :: automsk = .false.
    logical             :: wiener  = .false.
    logical             :: exists  = .false.
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! SETTERS
    procedure          :: reset_all
    procedure          :: reset_eos
    procedure, private :: reset_even
    procedure, private :: reset_odd
    procedure          :: reset_sum
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
    procedure          :: grid_fplane
    procedure          :: compress_exp
    procedure          :: sum_eos !< for merging even and odd into sum
    procedure          :: sum     !< for summing eo_recs obtained by parallel exec
    procedure          :: sampl_dens_correct_eos
    procedure          :: sampl_dens_correct_sum
    ! RECONSTRUCTION
    procedure          :: eorec
    ! DESTRUCTOR
    procedure          :: kill_exp
    procedure          :: kill
end type eo_reconstructor

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine new(self, p )
        use simple_fileio      , only: file_exists
        class(eo_reconstructor), intent(inout) :: self !< instance
        class(params), target,   intent(in)    :: p    !< parameters object (provides constants)
        type(image) :: imgtmp
        logical     :: neg  
        call self%kill
        ! set constants
        neg = .false.
        if( p%neg .eq. 'yes' ) neg = .true.
        self%box     = p%box
        self%smpd    = p%smpd
        self%nstates = p%nstates
        self%inner   = p%inner
        self%width   = p%width
        self%fny     = p%fny
        self%ext     = p%ext
        self%numlen  = p%numlen
        self%msk     = p%msk
        self%automsk = file_exists(p%mskfile)
        ! create composites
        if( self%automsk )then
            call self%envmask%new([p%box,p%box,p%box], p%smpd)
            call self%envmask%read(p%mskfile)
        endif
        call self%even%new([p%boxpd,p%boxpd,p%boxpd], p%smpd)
        call self%even%alloc_rho(p)
        call self%even%set_ft(.true.)
        call self%odd%new([p%boxpd,p%boxpd,p%boxpd], p%smpd)
        call self%odd%alloc_rho(p)
        call self%odd%set_ft(.true.)
        call self%eosum%new([p%boxpd,p%boxpd,p%boxpd], p%smpd)
        call self%eosum%alloc_rho(p, expand=.false.)
        ! set lfny
        call imgtmp%new([self%box,self%box,self%box], self%smpd)
        self%lfny = imgtmp%get_lfny(1)
        call imgtmp%kill
        ! set existence
        self%exists = .true.
    end subroutine new

    ! SETTERS

    !>  \brief  resets all
    subroutine reset_all( self )
        class(eo_reconstructor), intent(inout) :: self
        call self%reset_eos
        call self%reset_sum
    end subroutine reset_all

    !>  \brief  resets the even odd pairs
    subroutine reset_eos( self )
        class(eo_reconstructor), intent(inout) :: self
        call self%even%reset
        call self%odd%reset
    end subroutine reset_eos

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

    ! GETTERS

    !>  \brief  return the window functions used by eo_reconstructor
    function get_kbwin( self ) result( wf )
        use simple_kbinterpol, only: kbinterpol
        class(eo_reconstructor), intent(inout) :: self
        type(kbinterpol) :: wf
        wf = self%even%get_kbwin()
    end function get_kbwin

    !> \brief  for getting the resolution
    !> \param fsc05  target resolution a FSC=0.5
    !> \param fsc0143  target resolution a FSC=0.143
    subroutine get_res( self, fsc05, fsc0143 )
        class(eo_reconstructor), intent(in)  :: self !< instance
        real,                    intent(out) :: fsc05, fsc0143
        fsc0143 = self%fsc0143
        fsc05   = self%fsc05
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
        call self%even%write(trim(trim(adjustl(fbody))//'_even'//self%ext), del_if_exists=.true.)
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
        use simple_fileio      , only: file_exists
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
        use simple_fileio      , only: file_exists
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
    subroutine grid_fplane(self, o, fpl, pwght, mul, ran )
        use simple_ori, only: ori
        use simple_rnd, only: ran3
        class(eo_reconstructor), intent(inout) :: self  !< instance
        class(ori),              intent(inout) :: o     !< orientation
        class(image),            intent(inout) :: fpl   !< Fourier plane
        real, optional,          intent(in)    :: pwght !< external particle weight (affects both fplane and rho)
        real, optional,          intent(in)    :: mul   !< shift multiplication factor
        real, optional,          intent(in)    :: ran   !< external random number
        real    :: rran
        if( present(ran) )then
            rran = ran
        else
            rran = ran3()
        endif
        if( rran > 0.5 )then
            call self%even%inout_fplane(o, .true., fpl, pwght=pwght, mul=mul)
        else
            call self%odd%inout_fplane(o, .true., fpl, pwght=pwght, mul=mul)
        endif
    end subroutine grid_fplane

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

    !> \brief  for sampling density correction of the eo pairs
    subroutine sampl_dens_correct_eos( self, state, eonames )
        use simple_strings, only: int2str_pad
        use simple_fileio,  only: arr2file
        use simple_math,    only: get_resolution, calc_fourier_index
        use simple_masker,  only: masker
        class(eo_reconstructor), intent(inout) :: self       !< instance
        integer,                 intent(in)    :: state      !< state
        character(len=32),       intent(in)    :: eonames(2) !< even/odd filenames
        real, allocatable :: res(:), corrs(:)
        type(image)       :: even, odd
        type(masker)      :: volmasker
        integer           :: j, find
        ! make clipped volumes
        call even%new([self%box,self%box,self%box],self%smpd)
        call odd%new([self%box,self%box,self%box],self%smpd)
        ! correct for the uneven sampling density
        call self%even%sampl_dens_correct(maxits=1)
        call self%odd%sampl_dens_correct(maxits=1)
        ! reverse FT
        call self%even%bwd_ft
        call self%odd%bwd_ft
        ! clip
        call self%even%clip(even)
        call self%odd%clip(odd)
        if( self%automsk )then
            call even%mul(self%envmask)
            call odd%mul(self%envmask)
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
        ! write even/odd
        call even%write(trim(eonames(1)))
        call odd%write(trim(eonames(2)))
        ! forward FT
        call even%fwd_ft
        call odd%fwd_ft
        ! calculate FSC
        call even%fsc(odd, res, corrs)
        do j=1,size(res)
           write(*,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
        end do
        ! save, get & print resolution
        call arr2file(corrs, 'fsc_state'//int2str_pad(state,2)//'.bin')
        call get_resolution(corrs, res, self%fsc05, self%fsc0143)
        self%fsc05   = max(self%fsc05,self%fny)
        self%fsc0143 = max(self%fsc0143,self%fny)
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', self%fsc05
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', self%fsc0143
        ! the end
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
        call self%eosum%norm
        call self%eosum%clip(reference)
    end subroutine sampl_dens_correct_sum

    ! RECONSTRUCTION
    
    !> \brief  for reconstructing Fourier volumes according to the orientations 
    !!         and states in o, assumes that stack is open   
    subroutine eorec( self, fname, p, o, se, state, vol, mul, part, fbody )
        use simple_oris,            only: oris
        use simple_fileio,          only: file_exists
        use simple_sym,             only: sym
        use simple_params,          only: params
        use simple_gridding,        only: prep4cgrid
        use simple_imgfile,         only: imgfile, find_ldim_nptcls
        use simple_strings,         only: int2str_pad
        use simple_jiffys,          only: progress
        use simple_kbinterpol,      only: kbinterpol
        class(eo_reconstructor),    intent(inout) :: self      !< object
        character(len=*),           intent(in)    :: fname     !< spider/MRC stack filename
        class(params),              intent(in)    :: p         !< parameters
        class(oris),                intent(inout) :: o         !< orientations
        class(sym),                 intent(inout) :: se        !< symmetry element
        integer,                    intent(in)    :: state     !< state to reconstruct
        class(image),               intent(inout) :: vol       !< reconstructed volume
        real,             optional, intent(in)    :: mul       !< shift multiplication factor
        integer,          optional, intent(in)    :: part      !< partition (4 parallel rec)
        character(len=*), optional, intent(in)    :: fbody     !< body of output file
        type(image)       :: img, img_pad
        type(kbinterpol)  :: kbwin
        real              :: skewness
        integer           :: i, cnt, n, ldim(3), io_stat, filnum, state_glob
        integer           :: statecnt(p%nstates), alloc_stat, state_here
        character(len=32) :: eonames(2)
        call find_ldim_nptcls(fname, ldim, n)
        if( n /= o%get_noris() ) stop 'inconsistent nr entries; eorec; simple_eo_reconstructor'
        if( .not. present(part) )then
            if( p%eo .eq. 'aniso' ) stop 'eo=aniso not supported here, use simple_distr_exec!'
        endif
        kbwin = self%get_kbwin() 
        ! stash global state index
        state_glob = state
        ! make the images
        call img%new([p%box,p%box,1],p%smpd)
        call img_pad%new([p%boxpd,p%boxpd,1],p%smpd)
        ! even/odd partitioning
        if( o%get_nevenodd() == 0 ) call o%partition_eo('proj', [p%fromp,p%top])
        ! population balancing logics
        if( p%balance > 0 )then
            call o%balance( p%balance, NSPACE_BALANCE, p%nsym, p%eullims, skewness )
            write(*,'(A,F8.2)') '>>> PROJECTION DISTRIBUTION SKEWNESS(%):', 100. * skewness
        else
            call o%set_all2single('state_balance', 1.0)
        endif
        ! zero the Fourier volumes and rhos
        call self%reset_all
        write(*,'(A)') '>>> KAISER-BESSEL INTERPOLATION'
        statecnt = 0
        cnt      = 0
        do i=1,p%nptcls
            call progress(i, p%nptcls)
            if( i <= p%top .and. i >= p%fromp )then
                cnt = cnt+1
                state_here = nint(o%get(i,'state'))
                if( state_here > 0 .and. (state_here == state ) )then
                    statecnt(state) = statecnt(state)+1
                    call rec_dens
                endif
            endif
        end do
        ! undo fourier components expansion
        call self%compress_exp
        ! proceeds with density correction & output
        if( present(part) )then
            if( present(fbody) )then
                call self%write_eos(fbody//int2str_pad(state,2)//'_part'//int2str_pad(part,self%numlen))
            else
                call self%write_eos('recvol_state'//int2str_pad(state,2)//'_part'//int2str_pad(part,self%numlen))
            endif
        else
            if( present(fbody) )then
                eonames(1) = fbody//int2str_pad(state,2)//'_odd'//p%ext
                eonames(2) = fbody//int2str_pad(state,2)//'_even'//p%ext
                
            else
                eonames(1) = 'recvol_state'//int2str_pad(state,2)//'_odd'//p%ext
                eonames(2) = 'recvol_state'//int2str_pad(state,2)//'_even'//p%ext
            endif
            call self%sum_eos
            call self%sampl_dens_correct_eos(state, eonames)
            call self%sampl_dens_correct_sum(vol)
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
                use simple_rnd, only: ran3
                type(ori) :: o_sym, orientation
                integer   :: j, state, state_balance
                real      :: pw, eopart
                state         = nint(o%get(i, 'state'))
                state_balance = nint(o%get(i, 'state_balance'))
                if( state == 0 .or. state_balance == 0 ) return
                pw = 1.
                if( p%frac < 0.99 ) pw = o%get(i, 'w')
                if( pw > 0. )then
                    orientation = o%get_ori(i)
                    ! read image
                    if( p%l_distr_exec )then
                        call img%read(p%stk_part, cnt)
                    else
                        call img%read(fname, i)
                    endif
                    ! gridding
                    call prep4cgrid(img, img_pad, p%msk, kbwin)
                    ! e/o partitioning
                    eopart = ran3()
                    if( orientation%isthere('eo') )then
                        if( orientation%isevenodd() )eopart = orientation%get('eo')
                    endif
                    ! interpolation
                    if( p%pgrp == 'c1' )then
                        call self%grid_fplane(orientation, img_pad, pwght=pw, mul=mul, ran=eopart)
                    else
                        do j=1,se%get_nsym()
                            o_sym = se%apply(orientation, j)
                            call self%grid_fplane(o_sym, img_pad, pwght=pw, mul=mul)
                        end do
                    endif
                 endif
            end subroutine rec_dens
            
    end subroutine eorec

    ! DESTRUCTOR

    !>  \brief  is the expanded destructor
    subroutine kill_exp( self )
        class(eo_reconstructor), intent(inout)   :: self !< instance
        if( self%exists )then
            call self%even%dealloc_exp
            call self%odd%dealloc_exp
            call self%eosum%dealloc_exp
        endif
    end subroutine kill_exp

    !>  \brief  is a destructor
    subroutine kill( self )
        class(eo_reconstructor), intent(inout)   :: self !< instance
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
