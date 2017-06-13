program simple_test_parse
use simple_strings
use simple_defs
use simple_hash,  only: hash
use simple_chash, only: chash
use simple_sauron
implicit none
character(len=STDLEN) :: str
type(hash)            :: htab
type(chash)           :: chtab
str = 'micrograph=mymovie.mrc ctf_movie=mymovie_forctf.mrc e1=10.0 e2=20.0 e3=30.0 nparts=8 nthr=6'
call chtab%new(10)
call sauron_line_parser( str, htab, chtab )
call chtab%print_key_val_pairs
call htab%print
end program simple_test_parse
