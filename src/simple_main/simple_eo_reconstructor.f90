module simple_eo_reconstructor
use simple_reconstructor, only: reconstructor
use simple_image,         only: image
use simple_params,        only: params
use simple_defs           ! singleton
use simple_cmdline        ! singleton
use simple_jiffys         ! singleton
use simple_math           ! singleton
implicit none

public :: eo_reconstructor
private

type :: eo_reconstructor
    private
    type(reconstructor)    :: even
    type(reconstructor)    :: odd
    type(reconstructor)    :: eosum
    class(params), pointer :: pp=>null()
    character(len=4)       :: ext
    real                   :: fsc05, fsc0143, smpd, msk, fny, inner=0., width=10.
    integer                :: box=0, nstates=1, numlen=2, lfny=0
    logical                :: xfel=.false.
    logical                :: wiener=.false.
    logical                :: exists=.false.
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
    procedure          :: get_wfuns
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
    procedure          :: sum_eos ! for merging even and odd into sum
    procedure          :: sum     ! for summing eo_recs obtained by parallell exec
    procedure          :: sampl_dens_correct_eos
    procedure          :: sampl_dens_correct_sum
    procedure          :: eorec
    ! DESTRUCTOR
    procedure          :: kill
end type eo_reconstructor

contains

    ! CONSTRUCTOR
    
    !>  \brief  is a constructor
    subroutine new(self, p )
        use simple_ctf,    only: ctf
        class(eo_reconstructor),      intent(inout) :: self !< instance
        class(params), target,        intent(in)    :: p    !< parameters object (provides constants)
        type(image) :: imgtmp
        logical     :: neg  
        call self%kill
        ! set constants
        neg = .false.
        if( p%neg .eq. 'yes' ) neg = .true.
        self%pp      => p
        self%box     =  p%box
        self%smpd    =  p%smpd
        self%nstates =  p%nstates
        self%inner   =  p%inner
        self%width   =  p%width
        self%fny     =  p%fny
        self%ext     =  p%ext
        self%numlen  =  p%numlen
        self%msk     =  p%msk
        self%xfel    =  p%l_xfel
        ! create composites
        call self%even%new([p%boxpd,p%boxpd,p%boxpd], p%smpd, p%imgkind)
        call self%even%alloc_rho(p)
        call self%even%set_ft(.true.)
        call self%odd%new([p%boxpd,p%boxpd,p%boxpd], p%smpd, p%imgkind)
        call self%odd%alloc_rho(p)
        call self%odd%set_ft(.true.)
        call self%eosum%new([p%boxpd,p%boxpd,p%boxpd], p%smpd, p%imgkind)
        call self%eosum%alloc_rho(p)
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
    function get_wfuns( self ) result( wfs )
        use simple_winfuns, only: winfuns
        class(eo_reconstructor), intent(inout) :: self
        type(winfuns) :: wfs
        wfs = self%even%get_wfuns()
    end function get_wfuns
    
    !> \brief  for gettign the resolution
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
        character(len=*),        intent(in)    :: fbody
        call self%write_even(fbody)
        call self%write_odd(fbody)
    end subroutine write_eos
    
    !>  \brief  write the even reconstruction
    subroutine write_even( self, fbody )
        class(eo_reconstructor), intent(inout) :: self
        character(len=*),        intent(in)    :: fbody
        call self%even%write(trim(adjustl(fbody))//'_even'//self%ext, del_if_exists=.true.)
        call self%even%write_rho('rho_'//trim(adjustl(fbody))//'_even'//self%ext)
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
        inquire(FILE=even_vol, EXIST=here(1))
        inquire(FILE=even_rho, EXIST=here(2))
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
        inquire(FILE=odd_vol, EXIST=here(1))
        inquire(FILE=odd_rho, EXIST=here(2))
        if( all(here) )then
            call self%odd%read(odd_vol)
            call self%odd%read_rho(odd_rho)
        else
            call self%reset_odd
        endif
    end subroutine read_odd

    ! INTERPOLATION
    
    !> \brief  for gridding a Fourier plane
    subroutine grid_fplane( self, o, fpl, pwght, mul, ran, shellweights )
        use simple_ori, only: ori
        use simple_rnd, only: ran3
        class(eo_reconstructor), intent(inout) :: self            !< instance
        class(ori),              intent(inout) :: o               !< orientation
        class(image),            intent(inout) :: fpl             !< Fourier plane
        real, optional,          intent(in)    :: pwght           !< external particle weight (affects both fplane and rho)
        real, optional,          intent(in)    :: mul             !< shift multiplication factor
        real, optional,          intent(in)    :: ran             !< external random number
        real, optional,          intent(in)    :: shellweights(:) !< resolution weights
        real                                   :: rran
        if( present(ran) )then
            rran = ran
        else
            rran = ran3()
        endif
        if( rran > 0.5 )then
            call self%even%inout_fplane(o, .true., fpl, pwght=pwght, mul=mul, shellweights=shellweights)
        else
            call self%odd%inout_fplane(o, .true., fpl, pwght=pwght, mul=mul, shellweights=shellweights)
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
    
    !> \brief  for sampling density correction of the eo pairs
    subroutine sampl_dens_correct_eos( self, state )
        class(eo_reconstructor), intent(inout) :: self  !< instance
        integer,                 intent(in)    :: state !< state
        real, allocatable :: res(:), corrs(:)
        type(image)       :: even, odd
        integer           :: j
        ! make clipped volumes
        if( self%xfel )then
            call even%new([self%box,self%box,self%box],self%smpd,imgkind='xfel')
            call odd%new([self%box,self%box,self%box],self%smpd,imgkind='xfel')
        else
            call even%new([self%box,self%box,self%box],self%smpd)
            call odd%new([self%box,self%box,self%box],self%smpd)
        endif
        ! correct for the uneven sampling density
        call self%even%sampl_dens_correct
        call self%odd%sampl_dens_correct
        if( self%xfel )then
            ! no back transformation
        else
            ! reverse FT
            call self%even%bwd_ft
            call self%odd%bwd_ft
        endif
        ! clip
        call self%even%clip(even)
        call self%odd%clip(odd)
        if( self%xfel )then
            ! no masking or Fourier transformation
        else    
            if( self%inner > 1. )then
                call even%mask(self%msk, 'soft', inner=self%inner, width=self%width) 
                call odd%mask(self%msk, 'soft', inner=self%inner, width=self%width)
            else
                call even%mask(self%msk, 'soft') 
                call odd%mask(self%msk, 'soft')
            endif
            ! forward FT
            call even%fwd_ft
            call odd%fwd_ft
        endif
        ! calculate FSC
        call even%fsc(odd, res, corrs)
        call arr2file(corrs, 'fsc_state'//int2str_pad(state,2)//'.bin')
        do j=1,size(res)
           write(*,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
        end do
        ! get & print resolution
        call get_resolution(corrs, res, self%fsc05, self%fsc0143)
        self%fsc05   = max(self%fsc05,self%fny) 
        self%fsc0143 = max(self%fsc0143,self%fny)
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', self%fsc0143
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', self%fsc05
        ! the end
        deallocate(corrs, res)
        call even%kill
        call odd%kill
    end subroutine sampl_dens_correct_eos

    !> \brief  for sampling density correction, antialiasing, bwd_ft & normalization of the sum
    subroutine sampl_dens_correct_sum( self, reference )
        class(eo_reconstructor), intent(inout) :: self      !< instance
        class(image),            intent(inout) :: reference !< reference volume
        integer :: filtsz
        write(*,'(A)') '>>> SAMPLING DENSITY (RHO) CORRECTION & WIENER NORMALIZATION'
        call self%eosum%sampl_dens_correct
        if( self%xfel )then
            call self%eosum%clip(reference)
        else
            call self%eosum%bwd_ft
            call self%eosum%norm
            call self%eosum%clip(reference)
        endif
    end subroutine sampl_dens_correct_sum
    
    ! RECONSTRUCTION
    
    !> \brief  for reconstructing Fourier volumes according to the orientations 
    !!         and states in o, assumes that stack is open   
    subroutine eorec( self, fname, p, o, se, state, vol, mul, part, fbody, wmat )
        use simple_oris,            only: oris
        use simple_sym,             only: sym
        use simple_params,          only: params
        use simple_gridding,        only: prep4cgrid
        use simple_imgfile,         only: imgfile
        use simple_filterer,        only: resample_filter
        class(eo_reconstructor),    intent(inout) :: self      !< object
        character(len=*),           intent(in)    :: fname     !< spider/MRC stack filename
        class(params),              intent(in)    :: p         !< parameters
        class(oris),                intent(inout) :: o         !< orientations
        class(sym),                 intent(inout) :: se        !< symmetry element
        integer,                    intent(in)    :: state     !< state to reconstruct
        type(image),                intent(inout) :: vol       !< reconstructed volume
        real,             optional, intent(in)    :: mul       !< shift multiplication factor
        integer,          optional, intent(in)    :: part      !< partition (4 parallel rec)
        character(len=*), optional, intent(in)    :: fbody     !< body of output file
        real,             optional, intent(in)    :: wmat(:,:) !< shellweights
        type(image)       :: img, img_pad
        integer           :: i, cnt, n, ldim(3), filtsz, io_stat, filnum
        integer           :: filtsz_pad, statecnt(p%nstates), alloc_stat
        real, allocatable :: res(:), res_pad(:), wresamp(:)
        logical           :: doshellweight
        call find_ldim_nptcls(fname, ldim, n)
        if( n /= o%get_noris() ) stop 'inconsistent nr entries; eorec; simple_eo_reconstructor'
        doshellweight = present(wmat)
        ! make the images
        call img%new([p%box,p%box,1],p%smpd,p%imgkind)
        call img_pad%new([p%boxpd,p%boxpd,1],p%smpd,p%imgkind)
        ! make the stuff needed for shellweighting
        res        = img%get_res()
        res_pad    = img_pad%get_res()
        filtsz_pad = img_pad%get_filtsz()
        ! calculate weights
        if( p%frac < 0.99 ) call o%calc_hard_ptcl_weights(p%frac)
        ! zero the Fourier volumes and rhos
        call self%reset_all
        ! dig holds the state digit
        write(*,'(A)') '>>> KAISER-BESSEL INTERPOLATION'
        statecnt = 0
        call iterator(rec_dens)
        if( present(part) )then
            if( present(fbody) )then
                call self%write_eos(fbody//int2str_pad(state,2)//'_part'//int2str_pad(part,self%numlen))
            else
                call self%write_eos('recvol_state'//int2str_pad(state,2)//'_part'//int2str_pad(part,self%numlen))
            endif
        else
            call self%sum_eos
            call self%sampl_dens_correct_eos(state)
            call self%sampl_dens_correct_sum(vol)
        endif
        call img%kill
        call img_pad%kill
        ! report how many particles were used to reconstruct each state
        if( p%nstates > 1 )then
            write(*,'(a,1x,i3,1x,a,1x,i6)') '>>> NR OF PARTICLES INCLUDED IN STATE:', state, 'WAS:', statecnt(state)
        endif
        
        contains
        
            !> \brief  reconstruction iterator
            subroutine iterator( sub )
                interface
                    subroutine sub
                    end subroutine
                end interface
                cnt = 0
                do i=1,p%nptcls
                    call progress(i, p%nptcls)
                    if( i <= p%top .and. i >= p%fromp )then
                        cnt = cnt+1
                        if( nint(o%get(i,'state')) == state )then
                            statecnt(state) = statecnt(state)+1
                            call sub
                        endif
                    endif
                end do
            end subroutine iterator
        
            !> \brief  the density reconstruction functionality
            subroutine rec_dens
                use simple_ori, only: ori
                use simple_rnd, only: ran3
                type(ori) :: o_sym, orientation
                integer   :: j, state
                real      :: pw, ran
                state = nint(o%get(i, 'state'))
                if( state == 0 ) return
                pw = 1.
                if( p%frac < 0.99 ) pw = o%get(i, 'w')
                if( pw > 0. )then
                    orientation = o%get_ori(i)
                    if( p%l_distr_exec )then
                        call img%read(p%stk_part, cnt, p%l_xfel)
                    else
                        call img%read(fname, i, p%l_xfel)
                    endif
                    if( p%l_xfel )then
                        call img%pad(img_pad)
                    else
                        call prep4cgrid(img, img_pad, p%msk, wfuns=self%get_wfuns())
                    endif
                    if( doshellweight )then
                        wresamp = resample_filter(wmat(i,:), res, res_pad)
                    else
                        allocate(wresamp(filtsz_pad))
                        wresamp = 1.0
                    endif
                    if( p%pgrp == 'c1' )then
                        call self%grid_fplane(orientation, img_pad, pwght=pw, mul=mul, shellweights=wresamp)
                    else
                        ran = ran3()
                        do j=1,se%get_nsym()
                            o_sym = se%apply(orientation, j)
                            call self%grid_fplane(o_sym, img_pad, pwght=pw, mul=mul, ran=ran, shellweights=wresamp)
                        end do
                    endif
                    deallocate(wresamp)
                 endif
            end subroutine rec_dens
        
    end subroutine eorec
    
    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(eo_reconstructor), intent(inout)   :: self !< instance
        if( self%exists )then
            ! kill composites
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
