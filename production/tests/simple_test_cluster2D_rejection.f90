
program main
  use simple_defs,             only: COSMSKHALFWIDTH
  use simple_defs_fname,       only: MRC_EXT
  use simple_image,            only: image
  use simple_string,           only: string
  use simple_syslib,           only: simple_rename
  use simple_string_utils,     only: int2str
  use simple_sp_project,       only: sp_project
  use simple_fileio,           only: basename, swap_suffix, del_file
  use simple_imghead,          only: get_mrcfile_info
  use simple_stack_io,         only: stack_io
  use simple_imgarr_utils,     only: dealloc_imgarr, read_cavgs_into_imgarr
  use simple_cluster2D_rejector

  implicit none

  type :: test_sample
    type(string) :: projfile
    type(string) :: stkfile
    integer      :: msk
  end type test_sample
  
  type(image), allocatable :: cavg_imgs(:)
  integer,     allocatable :: states(:), cavg_inds(:)
  type(test_sample)        :: test_samples(15)
  type(sp_project)         :: spproj
  type(cluster2D_rejector) :: rejector
  type(string)             :: jpeg_out
  integer                  :: isample, xtiles, ytiles

  test_samples(1)%projfile = 'rejection_benchmarks/TMEM16A/microchunk_pass_1_1.simple'
  test_samples(1)%stkfile  = 'rejection_benchmarks/TMEM16A/microchunk_pass_1_1.mrc'
  test_samples(1)%msk      = 200

  test_samples(2)%projfile = 'rejection_benchmarks/TMEM16A/microchunk_pass_1_2.simple'
  test_samples(2)%stkfile  = 'rejection_benchmarks/TMEM16A/microchunk_pass_1_2.mrc'
  test_samples(2)%msk      = 200

  test_samples(3)%projfile = 'rejection_benchmarks/Betagal/microchunk_pass_1_1.simple'
  test_samples(3)%stkfile  = 'rejection_benchmarks/Betagal/microchunk_pass_1_1.mrc'
  test_samples(3)%msk      = 200

  test_samples(4)%projfile = 'rejection_benchmarks/Betagal/microchunk_pass_1_2.simple'
  test_samples(4)%stkfile  = 'rejection_benchmarks/Betagal/microchunk_pass_1_2.mrc'
  test_samples(4)%msk      = 200

  test_samples(5)%projfile = 'rejection_benchmarks/Betagal/microchunk_pass_1_3.simple'
  test_samples(5)%stkfile  = 'rejection_benchmarks/Betagal/microchunk_pass_1_3.mrc'
  test_samples(5)%msk      = 200

  test_samples(6)%projfile = 'rejection_benchmarks/PepT2/microchunk_pass_1_1.simple'
  test_samples(6)%stkfile  = 'rejection_benchmarks/PepT2/microchunk_pass_1_1.mrc'
  test_samples(6)%msk      = 160

  test_samples(7)%projfile = 'rejection_benchmarks/PepT2/microchunk_pass_1_2.simple'
  test_samples(7)%stkfile  = 'rejection_benchmarks/PepT2/microchunk_pass_1_2.mrc'
  test_samples(7)%msk      = 160

  test_samples(8)%projfile = 'rejection_benchmarks/PepT2/microchunk_pass_1_3.simple'
  test_samples(8)%stkfile  = 'rejection_benchmarks/PepT2/microchunk_pass_1_3.mrc'
  test_samples(8)%msk      = 160

  test_samples(9)%projfile = 'rejection_benchmarks/SLC/microchunk_pass_1_1.simple'
  test_samples(9)%stkfile  = 'rejection_benchmarks/SLC/microchunk_pass_1_1.mrc'
  test_samples(9)%msk      = 160

  test_samples(10)%projfile = 'rejection_benchmarks/SLC/microchunk_pass_1_2.simple'
  test_samples(10)%stkfile  = 'rejection_benchmarks/SLC/microchunk_pass_1_2.mrc'
  test_samples(10)%msk      = 160

  test_samples(11)%projfile = 'rejection_benchmarks/SLC/microchunk_pass_1_3.simple'
  test_samples(11)%stkfile  = 'rejection_benchmarks/SLC/microchunk_pass_1_3.mrc'
  test_samples(11)%msk      = 160

  test_samples(12)%projfile = 'rejection_benchmarks/Not1/microchunk_pass_1_1.simple'
  test_samples(12)%stkfile  = 'rejection_benchmarks/Not1/microchunk_pass_1_1.mrc'
  test_samples(12)%msk      = 160

  test_samples(13)%projfile = 'rejection_benchmarks/MotAB/microchunk_pass_1_1.simple'
  test_samples(13)%stkfile  = 'rejection_benchmarks/MotAB/microchunk_pass_1_1.mrc'
  test_samples(13)%msk      = 250

  test_samples(14)%projfile = 'rejection_benchmarks/CLC/microchunk_pass_1_1.simple'
  test_samples(14)%stkfile  = 'rejection_benchmarks/CLC/microchunk_pass_1_1.mrc'
  test_samples(14)%msk      = 209

  test_samples(15)%projfile = 'rejection_benchmarks/CLC/microchunk_pass_1_2.simple'
  test_samples(15)%stkfile  = 'rejection_benchmarks/CLC/microchunk_pass_1_2.mrc'
  test_samples(15)%msk      = 209

  do isample=1, size(test_samples)
    write(*,*) test_samples(isample)%projfile%to_char()
    call spproj%read(test_samples(isample)%projfile)
    call spproj%add_cavgs2os_out(test_samples(isample)%stkfile, spproj%get_smpd())
    cavg_imgs = read_cavgs_into_imgarr(spproj)

    call rejector%new(cavg_imgs, real(test_samples(isample)%msk))
    call rejector%reject_pop(spproj%os_cls2D)
    states = rejector%get_states()
    call spproj%os_cls2D%set_all('state', states)
    jpeg_out = swap_suffix(test_samples(isample)%stkfile, string('_pop_rejected.jpg'), string(MRC_EXT))
    call spproj%cavgs2jpg(cavg_inds, jpeg_out, xtiles, ytiles, ignore_states=.false., invert_states=.true.)

    call rejector%new(cavg_imgs, real(test_samples(isample)%msk))
    call rejector%reject_res(spproj%os_cls2D)
    states = rejector%get_states()
    call spproj%os_cls2D%set_all('state', states)
    jpeg_out = swap_suffix(test_samples(isample)%stkfile, string('_res_rejected.jpg'), string(MRC_EXT))
    call spproj%cavgs2jpg(cavg_inds, jpeg_out, xtiles, ytiles, ignore_states=.false., invert_states=.true.)

    call rejector%new(cavg_imgs, real(test_samples(isample)%msk))
    call rejector%reject_mask()
    states = rejector%get_states()
    call spproj%os_cls2D%set_all('state', states)
    jpeg_out = swap_suffix(test_samples(isample)%stkfile, string('_mask_rejected.jpg'), string(MRC_EXT))
    call spproj%cavgs2jpg(cavg_inds, jpeg_out, xtiles, ytiles, ignore_states=.false., invert_states=.true.)
    
    call rejector%new(cavg_imgs, real(test_samples(isample)%msk))
    call rejector%reject_local_variance()
    states = rejector%get_states()
    call spproj%os_cls2D%set_all('state', states)
    jpeg_out = swap_suffix(test_samples(isample)%stkfile, string('_local_variance_rejected.jpg'), string(MRC_EXT))
    call spproj%cavgs2jpg(cavg_inds, jpeg_out, xtiles, ytiles, ignore_states=.false., invert_states=.true.)

    cavg_imgs = read_cavgs_into_imgarr(spproj)
    call rejector%new(cavg_imgs, real(test_samples(isample)%msk))
    call rejector%reject_pop(spproj%os_cls2D)
    call rejector%reject_res(spproj%os_cls2D)
    call rejector%reject_mask()
    call rejector%reject_local_variance()
    states = rejector%get_states()
    call spproj%os_cls2D%set_all('state', states)
    jpeg_out = swap_suffix(test_samples(isample)%stkfile, string('_selected.jpg'), string(MRC_EXT))
    call spproj%cavgs2jpg(cavg_inds, jpeg_out, xtiles, ytiles, ignore_states=.false., invert_states=.false.)

    call rejector%kill()
    call spproj%kill()
    call dealloc_imgarr(cavg_imgs)

  end do

end program main
