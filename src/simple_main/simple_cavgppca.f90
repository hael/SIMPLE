!==Class simple_
!
module simple_cavgppca
use simple_defs              ! singleton
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_jiffys,           only: alloc_err
use simple_image,            only: image
use simple_build,            only: build
use simple_params,           only: params
use simple_ran_tabu,         only: ran_tabu
use simple_ppca,             only: ppca
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_matchpursuit,     only: matchpursuit

implicit none

public :: cavgppca
private

type cavgppca

  private
    class(params), pointer   :: pp=>null()
    class(build), pointer    :: pb=>null()
    type(ppca)               :: pca              ! Probabilistic PCA object
    type(oris), pointer      :: oris             ! Local oris (current class)
    type(image), pointer     :: ptr_ptcls(:)     ! Pointer to Particles stack
    real, allocatable        :: D(:,:)           ! Data matrix (ptcls)
    real, allocatable        :: avg(:)           ! average
    real, allocatable        :: evals(:)         ! Eigenvalues from PPCA
    real                     :: msk       = 0.
    integer                  :: n         = 0        ! number of ptcls
    integer                  :: nD        = 0        ! number of pixels/ptcl
    integer                  :: nfeats    = 10       ! number of latent vectors
    integer                  :: nselfeats = 10       ! number of selected latent vectors
    character(len=STDLEN)    :: metric    = 'realcc' ! Pixel of PFT space
    logical                  :: usemsk    = .false.  ! Whether to mask ptcls
    logical                  :: exists    = .false.

  contains
    ! CONSTRUCTORS
    procedure          :: new
    ! DOERS
    procedure, private :: ptcls2mat
    procedure, private :: make_featimg
    procedure, private :: write_mat
    procedure          :: exec
    procedure          :: get_nlatentvecs
    procedure          :: get_latentvecs
    ! DESTRUCTOR
    procedure          :: kill

end type

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    subroutine new( self, p, b, os, stk, nfeats, metric, msk )
        class(cavgppca), intent(inout)       :: self
        class(params), intent(in),target     :: p
        class(build), intent(in),target      :: b
        class(oris), intent(in), target      :: os
        class(image), intent(in), target     :: stk(:)
        integer, intent(in)                  :: nfeats
        character(len=*),intent(in)          :: metric
        logical, intent(in), optional        :: msk    
        integer                              :: alloc_stat, ldim(3)
        call self%kill
        if( present(msk) )self%usemsk=.true.
        self%pp         => p
        self%pb         => b
        self%oris       => os
        self%ptr_ptcls  => stk
        self%n          = os%get_noris()
        if( self%n<3 )stop 'Not enough data for ppca; in new;simple_cavg_ppca'
        self%metric    = metric
        self%nfeats    = nfeats
        self%nselfeats = nfeats
        write(*,'(A,I6)')'>>> PPCA SUBSPACE SIZE:',self%nfeats
        ! Data matrix
        if( self%usemsk )then
            self%msk = p%msk
            self%nD = self%ptr_ptcls(1)%get_npix( self%msk )
        else
            ldim    = self%ptr_ptcls(1)%get_ldim()
            self%nD = ldim(1)*ldim(2)*ldim(3)
        endif
        allocate(self%D(self%n,self%nD), self%evals(self%nfeats), &
            stat=alloc_stat)
        call alloc_err('In: new; simple_cavgppca 1', alloc_stat)
        ! average
        allocate(self%avg( self%nD ), stat=alloc_stat)
        call alloc_err('In: new; simple_cavgppca 2', alloc_stat)
        self%avg    = 0.
        self%exists = .true.
    end subroutine new

    subroutine write_mat( self )
        use simple_jiffys,   only: alloc_err, progress, get_fileunit
        class(cavgppca), intent(inout)   :: self
        integer :: i,j,filnum
        character(len=STDLEN) :: str
        filnum = get_fileunit()
        open(unit=filnum, status='REPLACE', action='WRITE', file='mat.txt')
        do i=1,self%n
            str = ''
            do j=1,self%n-1
                write(filnum,'(F12.6)',advance='no') self%D(i,j)
            enddo
            write(filnum,'(F12.6)') self%D(i,self%n)
        enddo
        close(filnum)
    end subroutine write_mat

    !>  \brief  is the master routine. Performs PPCA and determines 
    !>  the number of features based on plateauing of cumulative correlation
    subroutine exec( self )
        class(cavgppca), intent(inout)   :: self
        type(matchpursuit)               :: mpursuit
        integer :: featlim
        integer :: icls       ! debug
        real,allocatable :: featavg(:) ! debug
        character(len=STDLEN) :: fname ! debug
        ! extracts data & convert format
        call self%ptcls2mat
        !call self%write_mat
        ! PPCA
        call self%pca%new( self%n,self%nD,self%nfeats,Dmat=self%D,dorot=.true.,doprint=.false. )
        call self%pca%master( 'dummy1', 0, 'dummy2', 200 )
        self%evals     = self%pca%get_evals()
        ! Feature selection
        if( self%nfeats>1 )then
            ! matching pursuit
            call mpursuit%new( self%pca%get_loadings_ptr(), self%pca%get_feats_ptr(), self%avg, D=self%D, doprint=.true. )
            call mpursuit%exec( [1,self%nfeats], self%nselfeats )
            call mpursuit%kill
        else
            self%nselfeats = 1
        endif
        ! debug
        !self%nselfeats = self%nfeats
        ! end debug
        call self%make_featimg()
    end subroutine exec

    !> Backprojects features to images and alters orinal images stack
    subroutine make_featimg( self )
        class(cavgppca), intent(inout)   :: self
        real, allocatable                :: tmpvec(:)
        integer                          :: i, alloc_stat
        allocate( tmpvec(self%nD),stat=alloc_stat )
        call alloc_err('In: make_featimg; simple_cavgppca', alloc_stat)
        write(*,'(A,I3,A)')'>>> RESAMPLING USING ',self%nselfeats,' features'
        do i=1,self%n
            tmpvec = self%pca%generate( i,self%avg )
            if( self%usemsk )then
                self%ptr_ptcls(i) = 0.
                call self%ptr_ptcls(i)%serialize(tmpvec, self%msk)
                !call self%pb%img%mask( self%pp%msk,'hard' )
            else
                call self%ptr_ptcls(i)%serialize( pcavec=tmpvec )
            endif
        enddo
        deallocate( tmpvec )
    end subroutine make_featimg

    !> Formats image data into array
    subroutine ptcls2mat( self )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_stat,   only: normalize
        class(cavgppca), intent(inout) :: self
        real, allocatable              :: tmpD(:)
        integer                        :: i
        logical                        :: err
        ! normalization
        self%avg = 0.
        do i=1,self%n
            if( self%ptr_ptcls(i)%is_ft() )call self%ptr_ptcls(i)%bwd_ft
        enddo
        !$omp parallel do default(shared) private(i,err,tmpD)
        do i=1,self%n
            if( self%usemsk )then
                call self%ptr_ptcls(i)%serialize(tmpD, self%msk)
            else
                call self%ptr_ptcls(i)%serialize( pcavec=tmpD )
            endif
            err = .false.
            call normalize( tmpD,err )
            if( err )write(*,'(a,i7)') 'WARNING: variance zero! image nr: ', i
            self%D(i,:) = tmpD
            deallocate(tmpD)
        enddo
        !$omp end parallel do
        ! average
        do i=1,self%n
            self%avg = self%avg+self%D(i,:)/real(self%n)
        enddo
        ! centering
        do i=1,self%n
            self%D(i,:) = self%D(i,:)-self%avg
        enddo
    end subroutine ptcls2mat                

    ! GETTERS

    !>  \brief  retruns the number of selected features (latent vectors)
    function get_nlatentvecs( self )result( n )
        class(cavgppca), intent(inout) :: self
        integer  :: n
        n = self%nselfeats
    end function get_nlatentvecs

    !>  \brief  retruns the number of selected features (latent vectors)
    subroutine get_latentvecs( self, vecs )
        class(cavgppca), intent(inout)   :: self
        real, intent(inout), allocatable :: vecs(:,:)
        real, allocatable :: vec(:)
        integer           :: i,alloc_stat
        allocate( vecs(self%n,self%nselfeats),stat=alloc_stat )
        call alloc_err('In: new; simple_cavgppca 1', alloc_stat)
        do i=1,self%n
            vec       = self%pca%get_feat( i )
            vecs(i,:) = vec( 1:self%nselfeats )
            deallocate( vec )
        enddo
    end subroutine get_latentvecs

    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(cavgppca), intent(inout)   :: self
        if( self%exists )then
            self%pp        => null()
            self%oris      => null()
            self%ptr_ptcls => null()
            if( allocated(self%D) )deallocate(self%D)
            if( allocated(self%avg) )deallocate(self%avg)
            call self%pca%kill
            self%exists = .false.
        endif
    end subroutine kill
    
end module simple_cavgppca
