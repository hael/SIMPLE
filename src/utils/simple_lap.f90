! linear assignment problem solver, adapted from https://github.com/faolane/LAP/blob/master/lap.f90
module simple_lap
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
implicit none

public :: lap

private
type :: lap
    ! norm for matrix representation C(j,i) : j = columns, i = rows
    real,    allocatable :: CC(:,:)   ! cost matrix (min sum)
    integer, allocatable :: M(:,:)    ! mask matrix
    integer, allocatable :: rowCover(:), colCover(:) !cover row and cols
    integer :: n      ! dimension of C - assumed to be a (nxn) square matrix
    integer :: mode   ! 0: minimize sum, 1: maximize sum
    integer :: pathRow0, pathCol0  ! starting point for path finding part
contains
! CONSTRUCTOR
    procedure          :: new
    procedure          :: solve_lap
    procedure, private :: step1
    procedure, private :: step2
    procedure, private :: step3
    procedure, private :: step4
    procedure, private :: step5
    procedure, private :: step6
end type lap

contains
    ! CONSTRUCTORS
    subroutine new( self, mat )
        class(lap), target, intent(inout) :: self
        real,               intent(in)    :: mat(:,:)
        self%pathRow0 = 0
        self%pathCol0 = 0
        self%n        = size(mat, 1)
        allocate( self%CC(self%n,self%n) )
        self%CC = mat
    end subroutine new

    subroutine solve_lap(self, jSol)
        ! Implementation of the Munkres algorithm (also referred to as the Hungarian
        ! algorithm). J. Munkres, Journal of the SIAM 5, 1 (1957)
        ! The following implementation is based on
        ! http://csclab.murraystate.edu/%7Ebob.pilgrim/445/munkres.html
        class(lap), target, intent(inout) :: self
        integer,            intent(inout) :: jSol(:) ! solution indices
        integer :: step, i, j
        logical :: done
        done = .false.
        step = 1
        allocate(self%M(self%n,self%n))     ! mask matrix - contains starred zeros
        allocate(self%rowCover(self%n))     ! to keep track of covered rows
        allocate(self%colCover(self%n))     ! to keep track of covered columns
        do i = 1, self%n
            self%M(:,i)      = 0
            self%rowCover(i) = 0
            self%colCover(i) = 0
        enddo
        do while(.not. done)
            select case(step)
                case(1)
                    call self%step1(step)
                case(2)
                    call self%step2(step)
                case(3)
                    call self%step3(step)
                case(4)
                    call self%step4(step)
                case(5)
                    call self%step5(step)
                case(6)
                    call self%step6(step)
                case default ! done
                    do i = 1, self%n
                        do j = 1, self%n
                            if (self%M(j,i) == 1) jSol(i) = j
                        enddo
                    enddo
                    done = .true.
            end select
        enddo
        deallocate(self%M)
        deallocate(self%rowCover)
        deallocate(self%colCover)
    end subroutine solve_lap
    
    ! row reduction : for each row find the smallest value and substract it from
    ! all elements of that row. Go to step 2.
    subroutine step1( self, step )
        class(lap), target, intent(inout) :: self
        integer,            intent(out)   :: step
        integer :: i, j
        real    :: minVal
        do i = 1, self%n
            minVal = self%CC(1,i)
            do j = 1, self%n
                if( self%CC(j,i) < minVal ) minVal = self%CC(j,i)
            enddo
            self%CC(:,i) = self%CC(:,i) - minVal
        enddo
        step = 2
    end subroutine step1
    
    ! Search for zeros.
    ! Find a zero (Z) in the matrix. If no zeros has been previously starred in
    ! its row and column then star Z. Go to step 3.
    subroutine step2( self, step )
        class(lap), target, intent(inout) :: self
        integer,            intent(out)   :: step
        integer :: i, j
        do i = 1,self% n
            do j = 1, self%n
                if( self%CC(j,i) == 0 .and. self%rowCover(i) == 0 .and. self%colCover(j) == 0 )then
                    self%M(j,i) = 1
                    self%rowCover(i) = 1
                    self%colCover(j) = 1
                endif
            enddo
        enddo
        ! uncovers
        do i = 1, self%n
            self%rowCover(i) = 0
            self%colCover(i) = 0
        enddo
        step = 3
    end subroutine step2
    
    ! cover each column containing a starred zero. If n column are covered
    ! the starred zero describe an optimal assignment and we are done otherwise
    ! go to step 4.
    subroutine step3( self, step )
        class(lap), target, intent(inout) :: self
        integer,            intent(out)   :: step
        integer :: colCount, i, j
        colCount = 0
        do i = 1, self%n
            do j = 1, self%n
                ! if starred and column is uncovered
                if( self%M(j,i) == 1 .and. self%colCover(j) == 0 )then
                    self%colCover(j) = 1
                    colCount         = colCount + 1
                endif
            enddo
        enddo
        if (colCount == self%n) then
            step = 0
        else
            step = 4
        endif
    end subroutine step3
    
    ! Find a uncovered zero and prime it. If there is no starred zero in the row
    ! go to step 5. Otherwise, cover the row and uncover the column containing
    ! the starred zero. Continue until no uncovered zeros is left. Go to step 6.
    subroutine step4( self, step )
        class(lap), target, intent(inout) :: self
        integer,            intent(out)   :: step
        logical :: done, starInRow
        integer :: i, j, row, col
        done = .false.
        do while (.not. done)
            ! find an uncovered zero
            row       = 0
            col       = 0
            starInRow = .false.
            loop1: do i = 1, self%n
                loop2: do j = 1, self%n
                    if( self%CC(j,i) == 0 .and. self%rowCover(i) == 0 .and. self%colCover(j) == 0 )then
                        row = i
                        col = j
                        exit loop1
                    endif
                enddo loop2
            enddo loop1
            if (row == 0) then !no zero uncoverred left
                done = .true.
                step = 6
            else
                self%M(col,row) = 2 !primed zero
                ! search if there is a starred zero in the same row
                do j = 1, self%n
                    if (self%M(j,row) == 1) then
                        starInRow = .true.
                        col = j
                    endif
                enddo
                if (starInRow) then ! if there is a starred zero in line
                    self%rowCover(row) = 1
                    self%colCover(col) = 0
                else ! if no starred zero in line
                    done = .true.
                    step = 5
                    self%pathRow0 = row
                    self%pathCol0 = col
                endif
            endif
        enddo
    end subroutine step4
    
    ! Augmenting path algorithm: construct a serie of alternating primed and
    ! starred zeros as follows. Let Z0 be the uncoverd primed zero found in
    ! step 4. Let Z1 be the starred zero in the column of Z0 (if any).
    ! Let Z2 be the primed zero in the row of Z1 (there will always be one).
    ! Continue until the series terminates at a primed zero that has no starred
    ! zero in its column. Then unstar each starred zeros of the series, star
    ! each primed zeros of the series, erase all primes and uncover every line
    ! and columns. Return to step 3.
    subroutine step5( self, step )
        class(lap), target, intent(inout) :: self
        integer,            intent(out)   :: step
        logical :: done
        integer :: i, j, row, col, pathCount, path(2*self%n+1, 2)
        pathCount         = 1
        path(pathCount,1) = self%pathRow0
        path(pathCount,2) = self%pathCol0
        done              = .false.
        do while (.not. done)
            ! search for a starred zero in column
            row = 0
            col = path(pathCount,2)
            do i = 1, self%n
                if( self%M(col,i) == 1 ) row = i
            enddo
            if (row /= 0) then ! update path
                pathCount = pathCount + 1
                path(pathCount,1) = row
                path(pathCount,2) = path(pathCount-1,2)
            else
                done = .true.
            endif
            if( .not. done )then
                ! search for a prime zero in row
                do j = 1, self%n
                    if( self%M(j,row) == 2 ) col = j
                enddo
                ! update path
                pathCount = pathCount + 1
                path(pathCount,1) = path(pathCount-1,1)
                path(pathCount,2) = col
            endif
        enddo
        ! augment path
        do i = 1, pathCount
            if( self%M(path(i,2),path(i,1) ) == 1 )then
                self%M(path(i,2),path(i,1)) = 0
            else
                self%M(path(i,2),path(i,1)) = 1
            endif
        enddo
        ! clear covers and erase primes
        do i = 1, self%n
            self%rowCover(i) = 0
            self%colCover(i) = 0
            do j = 1, self%n
                if(self%M(j,i) == 2 ) self%M(j,i) = 0
            enddo
        enddo
        step = 3
    end subroutine step5
    
    ! Search for the smallest uncovered value and add it to covered rows
    ! and substract it from uncovered columns. Return to step 4.
    subroutine step6( self, step )
        class(lap), target, intent(inout) :: self
        integer,            intent(out)   :: step
        integer :: i, j
        real    :: minVal
        minVal = huge(minVal)
        do i = 1, self%n
            do j = 1, self%n
                if( self%rowCover(i) == 0 .and. self%colCover(j) == 0 .and. self%CC(j,i) < minVal )then
                    minVal = self%CC(j,i)
                endif
            enddo
        enddo
        do i = 1, self%n
            do j = 1, self%n
                if( self%rowCover(i) == 1 ) self%CC(j,i) = self%CC(j,i) + minVal
                if( self%colCover(j) == 0 ) self%CC(j,i) = self%CC(j,i) - minVal
            enddo
        enddo
        step = 4
    end subroutine step6
end module simple_lap