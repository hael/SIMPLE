module simple_align_pair
use simple_projector, only: projector
use simple_image,     only: image
use simple_polarft,   only: polarft
use simple_defs       ! singleton
implicit none

public :: align_pair, test_align_pair
private

type align_pair
    private
    real                        :: rot, trs, x, y, corr, msk
    integer                     :: sh, trsstep, irot, ring2
    integer, pointer            :: kfromto(:)=>null()
    type(image)                 :: img, tmp
    type(projector)             :: proj
    character(len=STDLEN)       :: space
    class(polarft), allocatable :: pimgs(:,:)
    class(polarft), pointer     :: refimg=>null()
    logical                     :: existence=.false.
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! GETTERS/SETTERS
    procedure :: exists
    procedure :: get_corr
    procedure :: set_refimg
    procedure :: set_pimg
    procedure :: reset_pimg
    procedure :: get_align_params
    procedure :: set_align_params
    procedure :: get_sh
    ! CALCULATORS
    procedure :: align
    procedure :: srch_all
    procedure :: srch_all_shifted
    procedure, private :: srch_rot
    procedure :: rocorr
    procedure, private :: srch_shifts
    procedure, private :: srch_refine_shift
    procedure, private :: srch_refine_shift_shifted
    procedure :: corr_ref_img
    ! DESTRUCTOR
    procedure :: kill
end type

contains

    !> \brief  is a constructor
    subroutine new( self, p, sh )
        use simple_params, only: params
        use simple_jiffys, only: alloc_err
        class(align_pair), intent(inout)  :: self
        class(params), intent(in), target :: p
        integer, intent(in), optional     :: sh
        integer                           :: alloc_stat
        if( self%existence ) call self%kill
        ! set constants
        self%trs   = p%trs
        self%msk   = p%msk
        self%ring2 = p%ring2
        self%sh    = nint(p%trs)
        if( present(sh) ) self%sh = sh
        self%trsstep = p%trsstep
        self%kfromto => p%kfromto
        ! make projector object
        self%proj = projector()
        if( self%sh == 0 )then
            allocate(self%pimgs(0:0,0:0), stat=alloc_stat)
        else
            allocate(self%pimgs(-self%sh:self%sh,-self%sh:self%sh), stat=alloc_stat)
        endif
        call alloc_err('new; simple_align_pair', alloc_stat)
        self%existence = .true.
    end subroutine
    
    ! GETTERS/SETTERS
    
    !> \brief  2 check existence of the object
    function exists( self ) result( existence )
        class(align_pair), intent(in) :: self
        logical :: existence
        existence = self%existence
    end function
    
    !> \brief  get correlation
    function get_corr( self ) result( corr )
        class(align_pair), intent(in) :: self
        real :: corr
        corr = self%corr
    end function
    
    !> \brief  set reference image
    subroutine set_refimg( self, refimg )
        class(align_pair), intent(inout)   :: self
        class(polarft), intent(in), target :: refimg
        if( refimg%exists() )then
            self%refimg => refimg
        else
            stop 'reference image does not exists; set_refimg; simple_align_pair'
        endif
    end subroutine
    
    !> \brief  set ptcl image
    subroutine set_pimg( self, img )
        class(align_pair), intent(inout), target :: self
        class(image), intent(in)                 :: img
        integer                                  :: ish, jsh
        call self%img%copy(img)
        if( self%sh > 0 )then
            ! create translation grid
            do ish=-self%sh,self%sh,self%trsstep
                do jsh=-self%sh,self%sh,self%trsstep
                    call self%tmp%copy(self%img)
                    call self%tmp%fwd_ft
                    call self%tmp%shift(real(ish),real(jsh))
                    call self%tmp%bwd_ft
                    call self%pimgs(ish,jsh)%new(self%kfromto, self%ring2)
                    call self%tmp%mask(self%msk,'soft')
                    call self%proj%img2polarft(self%tmp, self%pimgs(ish,jsh), self%msk)
                end do
            end do
        else
            call self%tmp%copy(self%img)
            call self%pimgs(0,0)%new(self%kfromto, self%ring2)
            call self%tmp%mask(self%msk,'soft')
            call self%proj%img2polarft(self%tmp, self%pimgs(0,0), self%msk) 
        endif
    end subroutine
    
    !> \brief  reset ptcl image (shifted or unshifted)
    subroutine reset_pimg( self, xy )
        class(align_pair), intent(inout), target :: self
        real, intent(in), optional               :: xy(2)
        if( present(xy) )then
            call self%tmp%copy(self%img)
            call self%tmp%fwd_ft
            call self%tmp%shift(xy(1),xy(2))
            call self%tmp%bwd_ft
            call self%pimgs(0,0)%new(self%kfromto, self%ring2)
            call self%tmp%mask(self%msk,'soft')   
            call self%proj%img2polarft(self%tmp, self%pimgs(0,0), self%msk)
        else
            call self%img%bwd_ft
            call self%tmp%copy(self%img)
            call self%pimgs(0,0)%new(self%kfromto, self%ring2)
            call self%tmp%mask(self%msk,'soft')
            call self%proj%img2polarft(self%tmp, self%pimgs(0,0), self%msk)
        endif
    end subroutine
    
    !> \brief  is for getting the alignment parameters
    function get_align_params( self ) result( rotxy )
        class(align_pair), intent(inout) :: self
        real                             :: rotxy(3)
        ! swap sign to fit convention
        rotxy(1) = 360.-self%rot
        if( rotxy(1) == 360. ) rotxy(1) = 0.
        rotxy(2) = -self%x
        rotxy(3) = -self%y
    end function
    
    !> \brief  is for setting the alignment parameters
    subroutine set_align_params( self, rotxy )
        class(align_pair), intent(inout) :: self
        real, intent(in) :: rotxy(3)
        ! swap sign to fit convention
        self%rot = 360.-rotxy(1)
        self%x   = -rotxy(2)
        self%y   = -rotxy(3)
    end subroutine
    
    !> \brief  is for getting the shift magnitude
    function get_sh( self ) result( sh )
        class(align_pair), intent(inout) :: self
        integer                          :: sh
        sh = self%sh
    end function
    
    ! CALCULATORS
    
    !> \brief  the high-level exhaustive alignment method
    function align( self ) result( rotxy )
        class(align_pair), intent(inout) :: self
        integer                          :: iang, ish, jsh
        real                             :: rotxy(3)
        call self%srch_all_shifted(iang, ish, jsh)
        if( self%sh > 0 )then
            call self%srch_refine_shift_shifted(ish, jsh)
        else
            self%x = 0.
            self%y = 0.
        endif
        rotxy = self%get_align_params() 
    end function

    !> \brief  to exhastively search the grid of translations and rotations
    subroutine srch_all( self, iang, ish_out, jsh_out, corr_out )
        class(align_pair), intent(inout) :: self
        integer, intent(out)             :: iang, ish_out, jsh_out
        real, intent(out), optional      :: corr_out
        real                             :: corr, rot
        integer                          :: ish, jsh, irot
        self%corr = -1.
        ! translation grid search
        do ish=-self%sh,self%sh,self%trsstep
            do jsh=-self%sh,self%sh,self%trsstep
                ! polar grid search
                call self%refimg%rotsrch(self%pimgs(ish,jsh), rot, corr, iang=irot)
                if( corr > self%corr )then  
                    self%corr = corr
                    self%rot  = rot
                    self%x    = real(ish)
                    self%y    = real(jsh)
                    self%irot = irot
                    iang      = irot
                    ish_out   = ish
                    jsh_out   = jsh
                endif
            end do
        end do
!         print *, 'shvec after grid srch: ', [self%x,self%y]
        if( present(corr_out) ) corr_out = self%corr
    end subroutine
    
    !> \brief  to exhastively search the grid of translations and rotations
    subroutine srch_all_shifted( self, iang, ish_out, jsh_out, corr_out )
        class(align_pair), intent(inout) :: self
        integer, intent(out)             :: iang, ish_out, jsh_out
        real, intent(out), optional      :: corr_out
        real                             :: corr, rot
        integer                          :: ish, jsh, irot
        self%corr = -1.
        ! translation grid search
        do ish=-self%sh,self%sh,self%trsstep
            do jsh=-self%sh,self%sh,self%trsstep
                ! polar grid search
                call self%refimg%shift([real(ish),real(jsh)])
                call self%refimg%rotsrch_shifted(self%pimgs(0,0), rot, corr, iang=irot)
                if( corr > self%corr )then  
                    self%corr = corr
                    self%rot  = rot
                    self%x    = real(ish)
                    self%y    = real(jsh)
                    self%irot = irot
                    iang      = irot
                    ish_out   = ish
                    jsh_out   = jsh
                endif
            end do
        end do
!         print *, 'shvec after grid srch: ', [self%x,self%y]
        if( present(corr_out) ) corr_out = self%corr
    end subroutine
    
    !> \brief  to exhastively search the grid of rotations
    subroutine srch_rot( self, ish_in, jsh_in, iang )
        class(align_pair), intent(inout) :: self
        integer, intent(in)              :: ish_in, jsh_in
        integer, intent(out)             :: iang
        real                             :: corr, rot
        integer                          :: irot
        call self%refimg%rotsrch(self%pimgs(ish_in,jsh_in), rot, corr, iang=irot)
        self%corr = corr
        self%rot  = rot
        self%irot = irot
        iang      = irot
    end subroutine

    !> \brief  to calculate a rotational correlation coefficient
    function rocorr( self, irot ) result( corr )
        class(align_pair), intent(inout) :: self
        integer, intent(in)              :: irot
        real                             :: corr
        if( self%refimg.eqdims.self%pimgs(0,0) )then
            corr = self%refimg%corr(self%pimgs(0,0), irot)
        else
            write(*,*) 'dim ref:', self%refimg%get_dims()
            write(*,*) 'dim ptcl:', self%pimgs(0,0)%get_dims()
            stop 'not equal dims; simple_align_pair; rocorr'
        endif
    end function
    
    !> \brief  to search translations only
    subroutine srch_shifts( self, iang )
        class(align_pair), intent(inout) :: self
        integer, intent(in)              :: iang
        real                             :: corr
        integer                          :: ish, jsh, ish_best, jsh_best
        self%corr = -1.
        ! translation grid search
        do ish=-self%sh,self%sh,self%trsstep
            do jsh=-self%sh,self%sh,self%trsstep
                corr = self%refimg%corr(self%pimgs(ish,jsh), iang)
                if( corr > self%corr )then  
                    self%corr = corr
                    self%x    = real(ish)
                    self%y    = real(jsh)
                    ish_best  = ish
                    jsh_best  = jsh
                endif
            end do
        end do
        call self%srch_refine_shift(ish_best, jsh_best)
    end subroutine
    
    !> \brief  refined shift search
    subroutine srch_refine_shift( self, ish_in, jsh_in )
        class(align_pair), intent(inout) :: self
        integer, intent(in)              :: ish_in, jsh_in
        real                             :: corr, xsh, ysh, lims(2,2)
        self%x = real(ish_in)
        self%y = real(jsh_in)
        lims(1,1) = max(-self%trs,self%x-0.75*real(self%trsstep))
        lims(1,2) = min( self%trs,self%x+0.75*real(self%trsstep))
        lims(2,1) = max(-self%trs,self%y-0.75*real(self%trsstep))
        lims(2,2) = min( self%trs,self%y+0.75*real(self%trsstep))
        self%corr = -1.
        xsh = lims(1,1)
        do while(xsh < lims(1,2))
            ysh = lims(2,1)
            do while(ysh < lims(2,2))
                ! shifted reference
                call self%reset_pimg([xsh,ysh])
                corr = self%refimg%corr(self%pimgs(0,0), self%irot)
                if( corr > self%corr )then
                    self%corr = corr
                    self%x    = xsh
                    self%y    = ysh
                endif
                ysh = ysh+0.25
            end do
            xsh = xsh+0.25
       end do
!         print *, 'shvec after refinement: ', [self%x,self%y]
    end subroutine
    
    !> \brief  refined shift search
    subroutine srch_refine_shift_shifted( self, ish_in, jsh_in )
        class(align_pair), intent(inout) :: self
        integer, intent(in)              :: ish_in, jsh_in
        real                             :: corr, xsh, ysh, lims(2,2)
        self%x = real(ish_in)
        self%y = real(jsh_in)
        lims(1,1) = max(-self%trs,self%x-0.75*real(self%trsstep))
        lims(1,2) = min( self%trs,self%x+0.75*real(self%trsstep))
        lims(2,1) = max(-self%trs,self%y-0.75*real(self%trsstep))
        lims(2,2) = min( self%trs,self%y+0.75*real(self%trsstep))
        self%corr = -1.
        xsh = lims(1,1)
        do while(xsh < lims(1,2))
            ysh = lims(2,1)
            do while(ysh < lims(2,2))
                ! shifted reference
                call self%refimg%shift([xsh,ysh])
                corr = self%refimg%corr_shifted(self%pimgs(0,0), self%irot)
                if( corr > self%corr )then
                    self%corr = corr
                    self%x    = xsh
                    self%y    = ysh
                endif
                ysh = ysh+0.25
            end do
            xsh = xsh+0.25
       end do
!         print *, 'shvec after refinement: ', [self%x,self%y]
    end subroutine
    
    !> \brief  calculate correlation to the 0,0 ptcl image given inplanes 
    function corr_ref_img( self, refinpl, ptclinpl ) result( r )
        class(align_pair), intent(inout) :: self
        integer, intent(in) :: refinpl, ptclinpl
        real :: r
        r = self%refimg%corr(self%pimgs(0,0), ptclinpl)
    end function
    
    ! DESTRUCTOR
    
    !> \brief  is a destructor
    subroutine kill( self )
        class(align_pair), intent(inout) :: self
        integer :: ish, jsh
        if( self%existence )then
            do ish=-self%sh,self%sh,self%trsstep
                do jsh=-self%sh,self%sh,self%trsstep
                    call self%pimgs(ish,jsh)%kill
                end do
            end do
            deallocate(self%pimgs)
            call self%img%kill
            call self%tmp%kill
            call self%refimg%kill
            self%existence = .false.
        endif
    end subroutine
    
    ! UNIT TEST
    
    !>  \brief  simple_align_pair unit test
    subroutine test_align_pair
        use simple_params, only: params
        use simple_rnd      ! singleton
        type(image)        :: img, img_rot, img_found, img_tmp
        type(align_pair)   :: ap 
        type(params)       :: p
        type(polarft)      :: pimg
        type(projector)    :: proj
        real               :: rotxy(3), corr, corrsum, rot, xsh, ysh, acorr, slask
        integer            :: i, ldim(3)
        logical            :: passed
        integer, parameter :: NROUNDS=1
        write(*,'(a)') '**info(simple_align_pair_unit_test): testing all functionality'
        passed = .false.
        call img%new([100,100,1],1.)
        call img%square(10)
        ldim = img%get_ldim()
        call img_rot%copy(img)
        call img_found%copy(img)
        proj = projector()
        corrsum = 0.
        call simulate
        acorr = corrsum/real(NROUNDS)
        if( acorr <= 0.98 )then
            print *, 'average correlation: ', acorr
            stop 'align_pair_test test failed'
        endif
        call img%kill
        call img_rot%kill
        call img_found%kill
        call ap%kill
        write(*,'(a)') 'SIMPLE_ALIGN_PAIR_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
        
        contains
        
            subroutine simulate
                do i=1,NROUNDS
                    rot = ran3()*360.
                    xsh = ran3()*5.-2.5
                    ysh = ran3()*5.-2.5
                    call img%rtsq(rot, xsh, ysh, img_rot)
                    p%trs = 3.
                    p%msk = real(ldim(1))/2.-2.            
                    p%kfromto = [2,20]
                    call ap%new(p)
                    call pimg%new(p%kfromto, nint(p%msk))
                    call proj%img2polarft(img, pimg, p%msk) ! THIS ONE BUGS OUT
                    call ap%set_refimg(pimg)
                    call ap%set_pimg(img_rot)
                    rotxy = ap%align()
                    call img_tmp%copy(img_rot)
                    call img_tmp%shift(-rotxy(2),-rotxy(3))
                    call img_tmp%rtsq(-rotxy(1), 0., 0., img_found)
                    ! calculate correlation
                    corr = img%corr(img_found, 20.)
                    corrsum = corrsum+corr
                end do
            end subroutine
            
    end subroutine
    
end module simple_align_pair
