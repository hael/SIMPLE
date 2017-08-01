! cgi_example.f90 --
!     Example of using CGI to implement an Internet service
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!
program cgi_example

    use cgi_protocol

    implicit none

    type(dict_struct), pointer :: dict => null() ! Initialisation is important!
    integer                    :: i
    integer                    :: luout
    integer                    :: steps
    real                       :: xmin
    real                       :: xmax
    real                       :: param
    real                       :: x, y
    character(len=20)          :: param_string

    !
    ! Get the CGI information
    ! and write the start of the HTML file (plus the start of the table)
    ! (Note: we include the text entry and submit button, so the user
    ! can easily select a new value)
    !
    call cgi_begin( output_html, dict, luout )

    call cgi_get( dict, "param", param )
    call cgi_get( dict, "param", param_string )

    write( luout, '(a)' ) '<html>'
    write( luout, '(a)' ) '<head><title>Table of function values</title></head>'
    write( luout, '(a)' ) '<body>'

    write( luout, '(a,a,a)' ) &
        'Parameter: <input type="text" name="param" value="', &
        trim(param_string), '"> <input type="submit">'

    write( luout, '(a)' ) '<table>'
    write( luout, '(a)' ) '<tr>'
    write( luout, '(3a)' ) '   <td>x</td><td>f(x) = exp(-x) * cos(ax)</td>'
    write( luout, '(a)' ) '</tr>'

    !
    ! Produce the table of function values
    !
    xmin  =   0.0
    xmax  =  10.0
    steps = 201

    do i = 1,steps
        x = (i-1) * 0.05

        y = func(param, x)

        write( luout, '(a)' ) '<tr>'
        write( luout, '(a,f10.4,a,f10.4,a)' ) &
            '    <td>', x, '</td><td>', y, '</td>'
        write( luout, '(a)' ) '</tr>'
    enddo

    write( luout, '(a)' ) '</table>'
    write( luout, '(a)' ) '</body>'
    write( luout, '(a)' ) '</html>'

    !
    ! We are done
    !
    call cgi_end

contains

real function func( param, x )
    real :: param, x

    func = exp(-x) * cos(param*x)

end function func

end program cgi_example
