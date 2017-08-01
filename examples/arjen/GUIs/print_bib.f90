! print_bib.f90 --
!     Demonstration program for reading bibliography XML files
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
program print_bib

    use xml_data_bibliography

    !
    ! Read in the entire file (leave out optional arguments)
    !
    call read_xml_file_bibliography( "example_bib.xml" )

    !
    ! Print the contents
    !
    do i = 1,size(book)
        write(*,'(a20,a,a)') book(i)%author, ' - ', trim(book(i)%title)
    enddo
end program print_bib
