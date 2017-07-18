!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Simple module: matching pursuit
!!  pearson correlation-based matching pursuit
module simple_matchpursuit
use simple_defs
use simple_jiffys, only: alloc_err

implicit none

public :: matchpursuit
private

type matchpursuit

  private
    real, pointer     :: pD(:,:) => null()     !< Pointer to centered data (NxD)
    real, pointer     :: pW(:,:) => null()     !< Pointer to loadings      (DxQ)
    real, pointer     :: pE(:,:) => null()     !< Pointer to expectations  (NxQ)
    real, pointer     :: pAVG(:) => null()     !< Pointer to data average  (D)
    integer                  :: N         = 0       !< number of observations
    integer                  :: Q         = 0       !< size of feature space
    integer                  :: D         = 0       !< dimension of observartions
    character(len=STDLEN)    :: fname               !< read from file
    logical                  :: doprint   = .false.
    logical                  :: inplace   = .true.
    logical                  :: exists    = .false.

  contains
    ! CONSTRUCTORS
    procedure          :: new
    ! DOERS
    procedure          :: exec
    procedure, private :: generate_cumul
    procedure, private :: get_corr
    ! DESTRUCTOR
    procedure          :: kill

end type matchpursuit

contains

    ! CONSTRUCTORS
    !>  \brief  is a constructor
    !>  should only be used if rotation onto orthogonal basis has been performed
    subroutine new( self, W, E, AVG, D, fname, doprint )
        class(matchpursuit), intent(inout) :: self 
        real, intent(in), pointer          :: W(:,:),E(:,:) !< loadings (DxQ) and expectations (NxQ)
        real, intent(in), target           :: AVG(:)     !< data average (D)
        real, intent(in), target, optional :: D(:,:)     !< arr for inplace matching (not to be used with fname)
        character( len=STDLEN),   optional :: fname      !< filename used for matching , disables inplace
        logical, intent(in),      optional :: doprint
        call self%kill
        self%N = size(D,1)
        self%D = size(D,2)
        self%Q = size(W,2)
        if( present(D) .and. .not.present(fname))then
            self%inplace = .true.
        else if( present(fname) .and. .not.present(D) )then
            self%inplace = .false.
            self%fname   = trim(fname)
        else
            stop 'Missing or inconstistent input arguments in simple_matchpursuit%new'
        endif
        if( present(doprint) )self%doprint=doprint
        if( self%inplace )then
            ! dimensions check
            if( self%N.ne.size(E,1) .or. self%D.ne.size(W,1) .or. self%Q.ne.size(E,2) .or. self%D.ne.size(AVG) ) &
                stop 'Inconstistent input matrix dimensions in simple_matchpursuit%new'
            ! pointer init
            self%pD     => D
        else
            ! stuff for file handle here ??
        endif
        self%pE     => E
        self%pW     => W
        self%pAVG   => AVG
        self%exists = .true.
    end subroutine new

    !>  \brief  performs pearson correlation-based matching pursuit & returns feature index
    subroutine exec( self, featlims, featlim, limit )
        class(matchpursuit), intent(inout)     :: self
        integer, intent(in)                    :: featlims(2) !< range of features to cycle trhough
        integer, intent(inout)                 :: featlim     !< feature index (return value)
        real,    intent(in), optional          :: limit       !< optional convergence criterion
        integer :: i, feat_ini
        real    :: avgcc, avgcc_prev
        real    :: lim
        if( self%doprint )write(*,'(A)')'>>> MATCHING PURSUIT'
        lim = 0.95
        if( present(limit) )lim=limit
        ! Init (featureless sampling if first is 1)
        feat_ini = featlims(1)-1
        call self%get_corr( [feat_ini,feat_ini], avgcc_prev )
        ! Iteration
        do featlim=featlims(1),featlims(2)
            call self%get_corr( [1,featlim], avgcc )
            if( avgcc_prev/avgcc>=lim )exit  ! convergence
            if( self%doprint )write(*,*)'Iteration',featlim,';cc=',avgcc,';cc_prev/cc=',avgcc_prev/avgcc
            avgcc_prev = avgcc
        enddo
        if( featlim>self%Q )featlim=self%Q
        write(*,'(A,F5.2,A,I3,A)')'>>> MATCHING PURSUIT CORRELATION=',avgcc,' (',featlim,' FEATURES)'
    end subroutine exec

    !>  \brief  calculates pearson correlation for a single feature range/value
    subroutine get_corr( self, feats, avgcc )
        use simple_stat,                    only: pearsn
        class(matchpursuit), intent(inout)     :: self
        real, intent(inout)      :: avgcc
        integer, intent(in)      :: feats(2)  !< resampling for subset features
        real    :: obs(self%D)    !< observed D-dimensional data
        real    :: model(self%D)  !< re-sampled
        real    :: X(self%D,1)
        integer :: i,file_stat,recsz,funit
        if( .not.self%inplace )then
            ! funit = get_fileunit()
            ! open(unit=self%funit, status='old', action='read', file=self%fname,&
            ! access='direct', form='unformatted', recl=recsz, iostat=file_stat)
            ! call fopen_err('master; simple_ppca, 1', file_stat)
        endif
        avgcc = 0.
        do i=1,self%N
            if( self%inplace )then
                obs   = self%pAVG+self%pD(i,:)   ! at each iteration to avoid cost of NxD matrix
            else
                ! read from file stuff here
                ! read(funit, rec=i) X(:,1) ! read from file
                ! obs = X(:,1)
            endif
            if( feats(1).eq.0 .and. feats(2).eq.0 )then
                model = self%pAVG  ! no resampling
            else
                model = self%generate_cumul(i,[feats(1),feats(2)])
            endif
            avgcc = avgcc+pearsn( obs,model )
        enddo
        avgcc = avgcc/real(self%N)
        !if( .not.self%inplace )close(unit=funit)
    end subroutine get_corr

    !>  \brief  is for sampling the generative model at a given observation index for a subset of features
    function generate_cumul( self, i, feats )result( dat ) 
        class(matchpursuit), intent(inout)  :: self
        integer, intent(in)        :: i, feats(2)
        real                       :: tmp(self%D,1),dat(self%D)
        dat = matmul(self%pW(:,feats(1):feats(2)),self%pE(i,feats(1):feats(2)))
        dat = dat+self%pAVG
    end function generate_cumul

    subroutine kill( self )
        class(matchpursuit), intent(inout)  :: self
        self%pD     => null()
        self%pE     => null()
        self%pW     => null()
        self%pAVG   => null()
        self%exists = .false.
    end subroutine kill

end module simple_matchpursuit
