module simple_roalgn
use simple_image,  only: image
use simple_oris,   only: oris
use simple_defs    ! singleton
implicit none
 
public :: roalgn, test_roalgn
private

type roalgn
    type(image)           :: img, img_rot
    class(oris), pointer  :: o_ptr=>null()
    integer, allocatable  :: pinds(:)
    character(len=STDLEN) :: stk
    real, allocatable     :: imgs(:,:), avg(:)
    real                  :: maxcorr=-1., msk
    integer               :: ldim(3)=[1,1,1], nptcls=0, npix=0
    logical               :: rnd=.false.
    logical               :: exists=.false.
  contains
    procedure          :: new
    procedure          :: align
    procedure, private :: read_ptcl
    procedure, private :: align_local
    procedure          :: kill
end type

interface roalgn
    module procedure constructor
end interface

integer, parameter :: MITS=7

contains

    !>  \brief  constructs the roalgn obj
    function constructor( stk, box, smpd, msk, o, rnd ) result( self )
        use simple_oris, only: oris
        type(roalgn)                    :: self
        character(len=*), intent(in)    :: stk
        integer, intent(in)             :: box
        real, intent(in)                :: smpd, msk
        class(oris), intent(in), target :: o
        character(len=*), intent(in)    :: rnd
        call self%new( stk, box, smpd, msk, o, rnd )
    end function

    !>  \brief  constructs the roalgn obj
    subroutine new( self, stk, box, smpd, msk, o, rnd )
        use simple_oris,   only: oris
        use simple_jiffys, only: alloc_err
        class(roalgn), intent(inout)           :: self
        character(len=*), intent(in)           :: stk
        integer, intent(in)                    :: box
        real, intent(in)                       :: smpd, msk
        class(oris), intent(in), target        :: o
        character(len=*), intent(in), optional :: rnd
        integer                                :: alloc_stat
        real, allocatable                      :: pcavec(:)
        call self%kill
        ! set constants:
        self%stk   = stk
        self%ldim  = [box,box,1]
        self%msk   = msk
        self%o_ptr => o
        call self%o_ptr%set_all('corr', -1.)
        if( present(rnd) )then
            if( rnd .eq. 'yes' )then
                self%rnd = .true.
            else 
                self%rnd = .false.
            endif
        endif
        self%nptcls = self%o_ptr%get_noris()
        ! make images
        call self%img%new([box,box,1],     smpd)
        call self%img_rot%new([box,box,1], smpd)
        call self%img%serialize(pcavec, self%msk)
        self%npix = size(pcavec)
        allocate( self%imgs(0:71,self%npix), self%avg(self%npix),&
        self%pinds(self%nptcls), stat=alloc_stat )
        call alloc_err('new; simple_roalgn', alloc_stat)
        self%exists = .true.       
    end subroutine
    
    !>  \brief  the master subroutine for the reference-free 2D alignment
    subroutine align( self )
        class(roalgn), intent(inout) :: self
        real    :: corr, corr_prev
        integer :: i, mode
        write(*,'(A)') '>>> REFERENCE-FREE ROTATIONAL ALIGNMENT'
        mode = 1
        call self%align_local(mode)
        do i=2,3
            corr_prev = self%o_ptr%get_avg('corr')
            mode = i
            call self%align_local(mode)
            corr = self%o_ptr%get_avg('corr')
            if( corr > corr_prev )then
                cycle 
            else
                exit
            endif
        end do
        ! swap the sign of the in-plane rotations (to fit convention)
        call self%o_ptr%e3swapsgn
    end subroutine
    
    !>  \brief  is for reading a particle image
    subroutine read_ptcl( self, i )
        class(roalgn), intent(inout) :: self
        integer, intent(in)          :: i
        real, allocatable            :: pcavec(:)
        integer                      :: r
        real                         :: theta
        call self%img%read(self%stk, i)
        do r=0,71
            theta = real(r)*5.
            call self%img%rtsq(theta, 0., 0., self%img_rot)
            call self%img_rot%serialize(pcavec, self%msk)
            self%imgs(r,:) = pcavec
            deallocate(pcavec)  
        end do
    end subroutine

    !>  \brief  this piece of magic is approximating the 3D in-plane angle
    subroutine align_local( self, mode )
        use simple_math,   only: hpsort, reverse
        use simple_jiffys, only: progress
        class(roalgn), intent(inout) :: self
        integer, intent(inout)       :: mode
        integer                      :: j, i, ii, irot1, irot2
        real                         :: corr, corr_prev, rot
        logical                      :: did_change
        call init
        write(*,'(A)') '>>> REFINEMENT'
        do j=1,MITS
            corr_prev = self%o_ptr%get_avg('corr')
            did_change = .false.
            call rnd_order
            do ii=1,self%nptcls
                i = self%pinds(ii)
                call progress(ii, self%nptcls)
                call self%read_ptcl(i)
                rot = self%o_ptr%e3get(i)
                irot1 = nint(rot/5.)             
                self%avg = self%avg-self%imgs(irot1,:)
                if( self%rnd )then
                    rot = rot_srch_rnd(i)
                else
                    rot = rot_srch()
                endif
                corr = self%o_ptr%get(i, 'corr')
                if( rot /= self%o_ptr%e3get(i) )then
                    call self%o_ptr%set(i, 'corr', self%maxcorr)
                    call self%o_ptr%e3set(i, rot)
                    irot2 = nint(rot/5.)
                    ! add back in new position
                    self%avg = self%avg+self%imgs(irot2,:)
                    did_change = .true.
                else
                    ! add back in old position
                    self%avg = self%avg+self%imgs(irot1,:)
                endif
            end do
            corr = self%o_ptr%get_avg('corr')
            write(*,"(1X,A,1X,I3,1X,A,1X,F7.4)") 'Iteration:', j, 'Correlation:', corr
            if( corr > corr_prev .and. did_change )then
                cycle 
            else
                exit
            endif
        end do

        contains
        
            function rot_srch( ) result( rot )
                use simple_stat, only: pearsn
                real    :: rot
                integer :: r
                self%maxcorr = -1.
                rot = 0.
                do r=0,71
                    corr = pearsn(self%avg, self%imgs(r,:))
                    if( corr > self%maxcorr )then
                        self%maxcorr = corr
                        rot = real(r)*5.
                    endif 
                end do
            end function
            
            function rot_srch_rnd( i ) result( rot )
                use simple_rnd, only: irnd_uni
                use simple_stat, only: pearsn
                integer, intent(in) :: i
                real    :: rot
                integer :: r, nbetter, indices(72), which
                real    :: corrs(72), corr
                do r=0,71
                    corrs(r+1) = pearsn(self%avg, self%imgs(r,:))
                    indices(r+1) = r
                end do
                call hpsort(72, corrs, indices)
                corr = self%o_ptr%get(i, 'corr')
                nbetter = 0
                do r=72,1,-1
                    if( corrs(r) >= corr )then
                        nbetter = nbetter+1
                        cycle
                    else
                        exit
                    endif
                end do
                if( nbetter > 1 )then
                    ! select a random better inplane
                    which = 72-irnd_uni(nbetter)+1 
                    rot = real(indices(which))*5.
                    self%maxcorr = corrs(which)
                else
                    ! select the best in-plane
                    rot = real(indices(72))*5.
                    self%maxcorr = corrs(72)
                endif               
            end function
            
            subroutine init
                integer           :: i, irot
                real              :: rot
                real, allocatable :: pcavec(:)
                write(*,'(A)') '>>> MAKING NEW 2D REFERENCE'
                call mk_pinds
                self%avg = 0.
                do i=1,self%nptcls
                    call progress(i, self%nptcls)
                    if( i == 1 )then
                        call self%img%read(self%stk, self%pinds(i))
                        call self%img%serialize(pcavec, self%msk)
                        self%avg = pcavec
                        deallocate(pcavec)
                    else
                        call self%read_ptcl(self%pinds(i))
                        rot = rot_srch()
                        call self%o_ptr%e3set(self%pinds(i), rot)
                        call self%o_ptr%set(self%pinds(i), 'corr', self%maxcorr)
                        irot = nint(rot/5.)
                        self%avg = self%avg+self%imgs(irot,:)
                    endif
                end do
                if( mode == 1 )then
                    mode = 2
                    call mk_pinds
                end if
            end subroutine

            subroutine mk_pinds
                real           :: corrs_copy(self%nptcls)
                integer        :: i
                if( mode == 1 )then ! random order
                    call rnd_order
                else                ! order according to corr
                    do i=1,self%nptcls
                        corrs_copy(i) = self%o_ptr%get(i,'corr')
                        self%pinds(i) = i
                    end do
                    call hpsort(self%nptcls, corrs_copy, self%pinds)
                    call reverse(self%pinds)
                endif
            end subroutine
            
            subroutine rnd_order
                use simple_ran_tabu, only: ran_tabu
                type(ran_tabu) :: rt
                integer        :: i
                rt = ran_tabu(self%nptcls)
                do i=1,self%nptcls
                    self%pinds(i) = rt%irnd()
                    call rt%insert(self%pinds(i))
                end do
                call rt%kill
            end subroutine

    end subroutine
    
    ! UNIT TEST
    
    subroutine test_roalgn
        use simple_image,  only: image
        use simple_oris,   only: oris
        use simple_rnd,    only: irnd_uni
        use simple_procimgfile, only: copy_imgfile, shrot_imgfile
        type(image)  :: img, img_rot
        type(roalgn) :: ralgn
        type(oris)   :: o
        integer      :: i, j
        real         :: corr
        logical      :: passed
        character(len=STDLEN), parameter :: fname='test_roalgn_squares.spi'
        character(len=STDLEN), parameter :: fname_aligned='test_roalgn_squares_aligned.spi'
        write(*,'(a)') '**info(simple_roalgn_unit_test): testing everything'
        passed = .false.
        call o%new(10)
        call img%new([100,100,1], 2.)
        call img_rot%new([100,100,1], 2.)
        call img%square(20)
        do i=1,10
            call img%rtsq(real(irnd_uni(72)-1)*5., 0., 0., img_rot)
            call img_rot%write(fname, i)
        end do
        call ralgn%new(fname, 100, 2., 45., o, 'yes')
        call ralgn%align
        call ralgn%kill
        call o%write('test_roalgn_algndoc.txt')
        call shrot_imgfile(fname, fname_aligned, o, 2.)
        corr = 0
        do i=1,9
            call img_rot%read(fname_aligned, i)
            do j=i+1,10   
                call img%read(fname_aligned, j)
                corr = corr+img%corr(img_rot)
            end do
        end do
        corr = corr/45.
        if( corr > 0.95 )then
            passed = .true.
        else
            print *, 'corr=', corr
        endif
        if( .not. passed )  stop 'simple_roalgn unit test failed'
        write(*,'(a)') 'SIMPLE_ROALGN_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(roalgn), intent(inout) :: self
        if( self%exists )then
            call self%img%kill
            call self%img_rot%kill
            deallocate(self%pinds, self%imgs, self%avg)
        endif
    end subroutine
    
end module simple_roalgn
