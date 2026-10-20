!@descr: filename helpers for refine2D class-average outputs
module simple_refine2D_fnames
use simple_defs_fname,   only: CAVGS_ITER_FBODY
use simple_string,       only: string
use simple_string_utils, only: int2str_pad, str2int
implicit none

private
public :: refine2D_cavgs_fbody
public :: refine2D_cavgs_fname
public :: refine2D_cavgs_half_fname
public :: refine2D_cavgs_ranked_fname
public :: refine2D_cavgs_iter

contains

    !> cavgs_iterNNN, the body of an iteration's class averages
    type(string) function refine2D_cavgs_fbody( iter ) result(fbody)
        integer, intent(in) :: iter
        fbody = string(CAVGS_ITER_FBODY)//int2str_pad(iter, 3)
    end function refine2D_cavgs_fbody

    !> the merged class averages of an iteration; ext is the extension with its dot (an image or JPG_EXT)
    type(string) function refine2D_cavgs_fname( iter, ext ) result(fname)
        integer,          intent(in) :: iter
        character(len=*), intent(in) :: ext
        fname = refine2D_cavgs_fbody(iter)
        fname = fname//ext
    end function refine2D_cavgs_fname

    !> the class averages of one half ('even' or 'odd') of an iteration
    type(string) function refine2D_cavgs_half_fname( iter, half, ext ) result(fname)
        integer,          intent(in) :: iter
        character(len=*), intent(in) :: half, ext
        fname = refine2D_cavgs_fbody(iter)
        fname = fname//'_'//trim(half)//ext
    end function refine2D_cavgs_half_fname

    !> the ranked class averages of an iteration
    type(string) function refine2D_cavgs_ranked_fname( iter, ext ) result(fname)
        integer,          intent(in) :: iter
        character(len=*), intent(in) :: ext
        fname = refine2D_cavgs_fbody(iter)
        fname = fname//'_ranked'//ext
    end function refine2D_cavgs_ranked_fname

    !> the iteration of a class-average file name: the digits from the end of the body to the extension; io_stat is
    !! nonzero when they do not form a number
    integer function refine2D_cavgs_iter( fname, ext, io_stat ) result(iter)
        class(string), intent(in)  :: fname, ext
        integer,       intent(out) :: io_stat
        integer :: first, last
        first = fname%substr_ind(CAVGS_ITER_FBODY) + len(CAVGS_ITER_FBODY)
        last  = fname%substr_ind(ext) - 1
        iter  = str2int(fname%to_char([first,last]), io_stat)
    end function refine2D_cavgs_iter

end module simple_refine2D_fnames
