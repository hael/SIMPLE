program simple_test_ansi_colors
  use simple_ansi_colors
  use simple_defs
  implicit none
  print '(a)', &
       color('Red',     c_red)     // NEWLINE // &
       color('Green',   c_green)   // NEWLINE // &
       color('Yellow',  c_yellow)  // NEWLINE // &
       color('Blue',    c_blue)    // NEWLINE // &
       color('Magenta', c_magenta) // NEWLINE // &
       color('Cyan',    c_cyan)    // NEWLINE // &
       color('White',   c_white)
end program simple_test_ansi_colors
