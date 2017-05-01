module simple_convergence_perptcl
use simple_oris,     only: oris
use simple_params,   only: params
use simple_cmdline,  only: cmdline
use simple_jiffys,   only: alloc_err
use simple_defs      ! use all in there
use simple_defs_conv ! use all in there
implicit none

public :: convergence_perptcl
private

real, parameter :: CONV_EPS = 0.5

type convergence_perptcl
    private
    class(oris),    pointer :: bap    => null() !< pointer to alignment oris object (a) part of build (b)
    class(params),  pointer :: pp     => null() !< pointer to parameters object
    class(cmdline), pointer :: pcline => null() !< pointer to command line object
    integer                 :: fromto(2)        !< particle range considered
    real, allocatable       :: mi_joint(:)      !< joint distribution overlaps
  contains
    ! updater (learning rate above)
    procedure          :: update_joint_distr_olap
    procedure          :: zero_joint_distr_olap
    ! generators
    procedure          :: gen_shift_larr
    procedure, private :: gen_larr
    procedure          :: identify_ptcls2process
    ! I/O
    procedure          :: write
    procedure          :: read
end type convergence_perptcl

interface convergence_perptcl
    module procedure constructor
end interface convergence_perptcl

contains

    function constructor( ba, p, cline ) result( self )
        class(oris),    target, intent(in) :: ba    !< alignment oris object (a) part of build (b)
        class(params),  target, intent(in) :: p     !< parameters object
        class(cmdline), target, intent(in) :: cline !< command line object
        type(convergence_perptcl) :: self
        integer :: alloc_stat
        self%bap       => ba
        self%pp        => p
        self%pcline    => cline
        self%fromto(1) =  p%fromp
        self%fromto(2) =  p%top
        allocate(self%mi_joint(p%fromp:p%top), stat=alloc_stat)
        call alloc_err("In: comple_convergence_perptcl :: constructor", alloc_stat)
        self%mi_joint  = 0.
    end function constructor

    subroutine update_joint_distr_olap( self )
        class(convergence_perptcl), intent(inout) :: self
        real, allocatable :: tmp(:)
        tmp = self%bap%get_all('mi_joint', self%fromto)
        self%mi_joint = (1. - CONV_EPS) * self%mi_joint + CONV_EPS * tmp
    end subroutine update_joint_distr_olap

    subroutine zero_joint_distr_olap( self )
        class(convergence_perptcl), intent(inout) :: self
        self%mi_joint = 0.
    end subroutine zero_joint_distr_olap

    function gen_shift_larr( self ) result( larr )
        class(convergence_perptcl), intent(in) :: self
        logical, allocatable :: larr(:)
        larr = self%gen_larr(FRAC_SH_LIM_PPTCL)
    end function gen_shift_larr

    function gen_larr( self, lim ) result( larr )
        class(convergence_perptcl), intent(in) :: self
        real,                       intent(in) :: lim
        logical, allocatable :: larr(:)
        integer :: alloc_stat
        allocate(larr(self%fromto(1):self%fromto(2)), stat=alloc_stat)
        call alloc_err("In: comple_convergence_perptcl :: gen_larr", alloc_stat)
        where( self%mi_joint >= lim )
            larr = .true.
        else where
            larr = .false.
        end where
    end function gen_larr

    function identify_ptcls2process( self ) result( inds )
        class(convergence_perptcl), intent(in) :: self
        logical, allocatable :: larr(:)
        integer, allocatable :: inds(:)
        integer :: ninds, iptcl, cnt, alloc_stat
        larr = self%gen_larr(MI_JOINT_LIM_2D)
        ninds = count(.not. larr)
        allocate( inds(ninds), stat=alloc_stat )
        call alloc_err("In: comple_convergence_perptcl :: identify_ptcls2process", alloc_stat)
        cnt = 0
        do iptcl=self%fromto(1),self%fromto(2)
            if( .not. larr(iptcl) )then
                cnt       = cnt + 1
                inds(cnt) = iptcl
            endif
        end do
    end function identify_ptcls2process

    subroutine write( self, fname )
        use simple_filehandling, only: get_fileunit
        class(convergence_perptcl), intent(in) :: self
        character(len=*),           intent(in) :: fname
        integer :: funit, io_stat
        funit = get_fileunit()
        open(unit=funit, status='REPLACE', action='WRITE', file=fname, access='STREAM')
        write(unit=funit,pos=1,iostat=io_stat) self%mi_joint
        ! Check if the write was successful
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(write): I/O error ',&
            io_stat, ' when writing to: ', trim(fname)
            stop 'I/O error; write; simple_convergence_perptcl'
        endif
        close(funit)
    end subroutine write

    subroutine read( self, fname )
        use simple_filehandling, only: get_fileunit
        class(convergence_perptcl), intent(inout) :: self
        character(len=*),           intent(in)    :: fname
        integer :: funit, io_stat
        funit = get_fileunit()
        open(unit=funit, status='OLD', action='READ', file=fname, access='STREAM')
        read(unit=funit,pos=1,iostat=io_stat) self%mi_joint
        ! Check if the read was successful
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(read): I/O error ',&
            io_stat, ' when reading from: ', trim(fname)
            stop 'I/O error; read; simple_convergence_perptcl'
        endif
        close(funit)
    end subroutine read

end module simple_convergence_perptcl
