# Where the small operations already live

Search here before writing a helper. Missing capability: extend the owner (an overload, an optional
argument, a stable variant), never a contained copy. A contained procedure is local policy or an adapter
over one of these (an empty-input guard, a fixed-length buffer), not a re-implementation.

| Need | Routine | Module |
| --- | --- | --- |
| ascending / descending index order of reals | `dists2order`, `scores2order` | `simple_srch_sort_loc` |
| descending order of DP values, ties kept in input order | `scores2order_stable` | `simple_srch_sort_loc` |
| sort values, or values with indices or paired weights | `hpsort` (generic) | `simple_srch_sort_loc` |
| reverse an array; k-th smallest; mask to indices | `reverse`, `selec`, `mask2inds` | `simple_srch_sort_loc` |
| first index of a value | intrinsic `findloc(arr, val, dim=1)` | |
| median (SP, DP); mean and sample deviation | `median`, `avg_sdev` | `simple_stat` |
| cross product (SP, DP, pure) | `cross` | `simple_math` |
| sub-sample peak offset from three samples | `parabolic_peak_offset` | `simple_math` |
| quote a shell argument | `posix_quote` | `simple_string_utils` |
| quote for Fortran/namelist text | `fortran_quote` | `simple_string_utils` |
| path from a root on, with a reusable position cache | `get_relative_path` (string, character) | `simple_fileio` |
| strip / read the last extension; body before a suffix | `strip_fext`, `fname2ext`, `get_fbody` | `simple_fileio` |
| file size by name (-1 absent); by open unit | `fname_size`, `funit_size` | `simple_fileio` |
| delete a file; remove a directory tree (links not followed) | `del_file`, `simple_rmdir` | `simple_syslib` |
| allocate / destroy an image array | `alloc_imgarr`, `dealloc_imgarr` | `simple_imgarr_utils` |
| the class-average CNN (scripts/sscnn): weights file, forward pass, group probabilities | `cavg_cnn` | `simple_cavg_cnn` |
| a class average and its halves as that CNN's input (numpy's Fourier resampling exactly, any size, then its scaling) | `cavg_cnn_input` | `simple_cavg_cnn_input` |
| does a window (corner coord, box) fit in the image | `window_inside` | `simple_image` |
| stop on a stack header with no images or a zero dimension | `check_stack_dims` | `simple_imghead` |
| elapsed seconds without touching the shared timer (threads) | `clock_count`, `clock_seconds_since` | `simple_timer` |
| one RGB pixel for `write_jpg(colorspec=3)` | `jpg_rgb_code` | `simple_jpg` |
| particles in a state; of those, updated at least once | `count_state_gt_zero`, `count_state_updated` | `oris` |

Same name is not same contract: check the convention (sample or population deviation, inclusive
bounds, tie order, precision, empty input) before replacing a local routine, and keep independent test
oracles independent of the routine they test.
