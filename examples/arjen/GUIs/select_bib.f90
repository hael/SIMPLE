! select_bib.f90 --
!     Demonstration program for selecting books by a single author from
!     bibliography XML files
!     Process the file as we go along
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
module bibliography_scan
    implicit none

    !
    ! Work arrays for storing the information from the XML file
    !
    character(len=20), dimension(2,100) :: attribs
    character(len=80), dimension(100)   :: data

    !
    ! Variables for storing the information that is to be printed
    !
    character(len=20) :: author
    logical           :: store = .false.

contains
subroutine startfunc( tag, attribs, error )
    character(len=*)                 :: tag
    character(len=*), dimension(:,:) :: attribs
    logical                          :: error

    ! Dummy - has no function in this case
end subroutine startfunc
subroutine endfunc( tag, error )
    character(len=*)                 :: tag
    logical                          :: error

    ! Dummy - has no function in this case
end subroutine endfunc

subroutine datafunc( tag, data, error )
    character(len=*)                 :: tag
    character(len=*), dimension(:)   :: data
    logical                          :: error

    if ( tag == "author" .and. index( data(1), "Swift") > 0 ) then
        author = data(1)
        store  = .true.
    endif
    if ( tag == "title" .and. store ) then
        write(*,'(a20,a,a)') author, ' - ', trim(data(1))
        store = .false.
    endif
end subroutine datafunc
end module bibliography_scan

program select_bib
    use xmlparse
    use bibliography_scan

    implicit none

    integer :: lunrep
    logical :: error

    !
    ! Read in the entire file (leave out optional arguments)
    !
    lunrep = 10
    open( lunrep, file = "select_bib.log" )
    call xml_process( "example_bib.xml", attribs, data, startfunc, &
             datafunc, endfunc, lunrep, error )

end program select_bib
