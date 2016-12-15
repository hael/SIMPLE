# plotgui.tcl --
#     Very simple GUI:
#     - Get one value
#     - Run the computational program
#     - Display the result
#
#     Example belonging to "Modern Fortran in Practice" by Arjen Markus
#
#     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
#     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
#     or send a letter to:
#     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
#

#
# Load a plotting package
#
package require Plotchart

#
# Create the user-interface elements
#
::ttk::frame  .toprow
::ttk::button .toprow.b -text Plot -command {putValue} -width 10
::ttk::label  .toprow.l -text Parameter:
::ttk::entry  .toprow.e -textvariable parameter

canvas .c -width 500 -height 400 -bg white

#
# Arrange them in the main window (.)
grid  .toprow.l .toprow.e .toprow.b -padx 5 -sticky w
grid  .toprow
grid  .c -  -

#
# Auxiliary procedures
#
proc putValue {} {
    global parameter
    global program
    global plot

    puts $program $parameter
    flush $program

    #
    # Clean up the graph
    #
    $plot plot data {} {}
    .c delete data
}

proc readData {channel} {
    global plot

    if { ![eof $channel] } {

        gets $channel line

        scan $line "%f %f" x y
        $plot plot data $x $y

    } else {
        close $channel
    }
}

#
# Set up the plot and start the program
#
set parameter 1.0

set plot [::Plotchart::createXYPlot .c {0.0 10.0 2.0} {-1.0 1.0 0.25}]

set program [open "|runprogram" r+]
fileevent $program readable [list readData $program]

#
# The event loop starts automatically ...
# We can just wait now
