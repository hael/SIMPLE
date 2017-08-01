module simple_orialgn_pair
use simple_ori,         only: ori
use simple_oris,        only: oris
use simple_simplex_opt, only: simplex_opt
use simple_opt_spec,    only: opt_spec
implicit none

public :: orialgn_pair, test_orialgn_pair
private

type orialgn_pair
    private
    type(ori)            :: dock_ori        !< the orientation relating ref with target
    class(oris), pointer :: reforis=>null() !< orientations of refernce
    type(oris)           :: targoris        !< orientations of target
    type(simplex_opt)    :: opt             !< optimizer
    type(opt_spec)       :: ospec           !< optimizer specification
    logical              :: mirror=.false.  !< to indicate mirror or not
    logical              :: exists=.false.  !< to indicate existence
  contains
    procedure :: new
    procedure :: get_docking
    procedure :: transform_oris
    procedure :: set_reforis
    procedure :: set_targoris
    procedure, private :: assign
    generic   :: assignment(=) => assign
    procedure :: dock
    procedure :: kill
end type

interface orialgn_pair
    module procedure :: constructor
end interface

real, parameter :: TOL=1e-10

contains
    
    !> \brief  is a constructor
    function constructor( reforis, targoris ) result( self )
        class(oris), intent(in), target, optional :: reforis
        class(oris), intent(in), optional         :: targoris
        type(orialgn_pair)                        :: self
        call self%new( reforis, targoris )
    end function

    !> \brief  is a constructor
    subroutine new( self, reforis, targoris )
        class(orialgn_pair), intent(inout)        :: self
        class(oris), intent(in), target, optional :: reforis
        class(oris), intent(in), optional         :: targoris
        real :: limits(3,2)
        limits = 0.
        limits(1,2) = 360.
        limits(2,2) = 180.
        limits(3,2) = 360.
        call self%kill
        call self%dock_ori%new
        if( present(reforis) )  call self%set_reforis(reforis)
        if( present(targoris) ) call self%set_targoris(targoris)
        call self%ospec%specify('simplex', 3, ftol=TOL,&
        limits=limits, cyclic=[.true.,.true.,.true.], nrestarts=5)
        call self%opt%new(self%ospec)
        self%exists=.false.
    end subroutine
    
    !> \brief  is a getter
    subroutine get_docking( self, o, mirr )
        class(orialgn_pair), intent(in) :: self
        class(ori), intent(out)         :: o
        logical, intent(out)            :: mirr
        o    = self%dock_ori
        mirr = self%mirror
    end subroutine
    
    !> \brief  transforms inoput oris according to docking
    subroutine transform_oris( self, os )
        class(orialgn_pair), intent(in) :: self
        class(oris), intent(inout)      :: os
        if( self%mirror ) call os%mirror3d
        call os%rot(self%dock_ori)
    end subroutine
    
    !> \brief  is a setter
    subroutine set_reforis( self, reforis )
        class(orialgn_pair), intent(inout) :: self
        class(oris), intent(in), target    :: reforis
        self%reforis => reforis
    end subroutine
    
    !> \brief  is a setter
    subroutine set_targoris( self, targoris )
        class(orialgn_pair), intent(inout) :: self
        class(oris), intent(in)            :: targoris
        self%targoris = targoris
    end subroutine
    
    !>  \brief is polymorphic assignment, overloaded as (=)
    subroutine assign( self_out, self_in )
        class(orialgn_pair), intent(inout) :: self_out
        class(orialgn_pair), intent(in)    :: self_in 
        call self_out%new
        self_out%reforis => self_in%reforis
        self_out%targoris = self_in%targoris
    end subroutine
    
    !> \brief  does the work
    subroutine dock( self, cost_out, iniori )
        use simple_rnd, only: ran3
        class(orialgn_pair), intent(inout), target :: self
        real, intent(out)                          :: cost_out
        class(ori), intent(inout), optional        :: iniori
        type(orialgn_pair) :: global_obj    
        type(ori)          :: global_ori   
        type(oris)         :: global_oris
        real               :: cost_lowest(2), dock_eul(2,3)
        integer            :: i
        logical            :: mirror
        call self%ospec%set_costfun(dock_cost)
        call global_ori%new
        do i=1,2
            mirror = .false.
            if( i == 2 ) mirror = .true.
            global_obj = self 
            call global_ori%new
            if( present(iniori) )then
                self%ospec%x = iniori%get_euler()
            else
                self%ospec%x(1) = ran3()*360.
                self%ospec%x(2) = ran3()*180.
                self%ospec%x(3) = ran3()*360.
            endif
            call self%opt%minimize(self%ospec, cost_lowest(i))
            dock_eul(i,:) = self%ospec%x
        end do
        if( cost_lowest(2) < cost_lowest(1) )then
            call self%dock_ori%set_euler(dock_eul(2,:))
            self%mirror = .true.
            cost_out = cost_lowest(2)
        else
            call self%dock_ori%set_euler(dock_eul(1,:))
            self%mirror = .false.
            cost_out = cost_lowest(1)
        endif
        call global_obj%kill
        call global_ori%kill
        call global_oris%kill
        
        contains
        
            !>  \brief  calculates the target vs. reference distance
            function dock_cost( vec, D ) result( cost )
                integer, intent(in) :: D
                real, intent(in)    :: vec(D)
                real                :: cost
                global_oris = global_obj%targoris
                if( mirror ) call global_oris%mirror3d
                call global_ori%set_euler([vec(1),vec(2),vec(3)])
                call global_oris%rot(global_ori)
                cost = global_obj%reforis.inpldist.global_oris
            end function

    end subroutine
    
    !> \brief  is the orialgn_pair unit test
    subroutine test_orialgn_pair
        use simple_rnd ! singleton
        type(oris)         :: o1, o2
        type(ori)          :: o
        type(orialgn_pair) :: oap
        integer            :: i
        character(len=20)  :: file
        character(len=2)   :: dig
        real               :: cost
        o1 = oris(20)
        call o1%spiral
        call o%new
        do i=1,10
            o2 = o1
            call o%set_euler([ran3()*360.,ran3()*180.,ran3()*360.])
            if( ran3() < 0.5 ) call o2%mirror3d
            call o2%rot(o) 
            oap = orialgn_pair(o1,o2)
            call oap%dock(cost)
            call oap%transform_oris(o2)
            write(dig,'(I2)') i
            file = 'oridoc_'//trim(adjustl(dig))//'.txt'
            call o2%write(file)
        end do
    end subroutine
    
    !> \brief  is a destructor
    subroutine kill( self )
        class(orialgn_pair), intent(inout) :: self
        call self%dock_ori%kill
        self%reforis => null()
        call self%targoris%kill
        call self%opt%kill
        call self%ospec%kill
    end subroutine
    
end module simple_orialgn_pair
