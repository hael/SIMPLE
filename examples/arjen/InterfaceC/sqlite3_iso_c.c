/* sqlite3_iso_c.c --

       Example belonging to "Modern Fortran in Practice" by Arjen Markus

       This work is licensed under the Creative Commons Attribution 3.0 Unported License.
       To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
       or send a letter to:
       Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

       Dummy C routine to illustrate the iso_c_binding module
       Accompanies test_iso_c.f90
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int sqlite3_column_count( void *stmt ) {

    return 3;
}

int Sqlite3DoC( void *handle, char *command, char *errmsg, int len_errmsg ) {

    printf( "Length: %d\n", strlen(command) );
    strcpy( errmsg, command );
    return len_errmsg;
}
