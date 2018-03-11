program simple_test_ansi_format_strs
  use simple_ansi_ctrls
  use simple_defs
  implicit none
  print '(a)', &
       format_str('Red',     C_RED)     // NEWLINE // &
       format_str('Green',   C_GREEN)   // NEWLINE // &
       format_str('Yellow',  C_YELLOW)  // NEWLINE // &
       format_str('Blue',    C_BLUE)    // NEWLINE // &
       format_str('Magenta', C_MAGENTa) // NEWLINE // &
       format_str('Cyan',    C_CYAN)    // NEWLINE // &
       format_str('White',   C_WHITE)
end program simple_test_ansi_format_strs
