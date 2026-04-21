
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
  type(test_sample)        :: test_samples(87)
  type(sp_project)         :: spproj
  type(cluster2D_rejector) :: rejector
  type(string)             :: jpeg_out
  integer                  :: isample, xtiles, ytiles


  ! TMEM16A
  test_samples(1)%projfile = 'rejection_benchmarks/TMEM16A/microchunk_pass_1_1.simple'
  test_samples(1)%stkfile  = 'rejection_benchmarks/TMEM16A/microchunk_pass_1_1.mrc'
  test_samples(1)%msk      = 200

  test_samples(2)%projfile = 'rejection_benchmarks/TMEM16A/microchunk_pass_1_2.simple'
  test_samples(2)%stkfile  = 'rejection_benchmarks/TMEM16A/microchunk_pass_1_2.mrc'
  test_samples(2)%msk      = 200

  test_samples(3)%projfile = 'rejection_benchmarks/TMEM16A/microchunk_pass_1_3.simple'
  test_samples(3)%stkfile  = 'rejection_benchmarks/TMEM16A/microchunk_pass_1_3.mrc'
  test_samples(3)%msk      = 200

  test_samples(4)%projfile = 'rejection_benchmarks/TMEM16A/microchunk_pass_2_1.simple'
  test_samples(4)%stkfile  = 'rejection_benchmarks/TMEM16A/microchunk_pass_2_1.mrc'
  test_samples(4)%msk      = 200

  test_samples(5)%projfile = 'rejection_benchmarks/TMEM16A/microchunk_pass_2_2.simple'
  test_samples(5)%stkfile  = 'rejection_benchmarks/TMEM16A/microchunk_pass_2_2.mrc'
  test_samples(5)%msk      = 200

  test_samples(6)%projfile = 'rejection_benchmarks/TMEM16A/microchunk_pass_2_3.simple'
  test_samples(6)%stkfile  = 'rejection_benchmarks/TMEM16A/microchunk_pass_2_3.mrc'
  test_samples(6)%msk      = 200

  test_samples(7)%projfile = 'rejection_benchmarks/TMEM16A/refchunk.simple'
  test_samples(7)%stkfile  = 'rejection_benchmarks/TMEM16A/refchunk.mrc'
  test_samples(7)%msk      = 200

  ! Betagal
  test_samples(8)%projfile = 'rejection_benchmarks/Betagal/microchunk_pass_1_1.simple'
  test_samples(8)%stkfile  = 'rejection_benchmarks/Betagal/microchunk_pass_1_1.mrc'
  test_samples(8)%msk      = 200

  test_samples(9)%projfile = 'rejection_benchmarks/Betagal/microchunk_pass_1_2.simple'
  test_samples(9)%stkfile  = 'rejection_benchmarks/Betagal/microchunk_pass_1_2.mrc'
  test_samples(9)%msk      = 200

  test_samples(10)%projfile = 'rejection_benchmarks/Betagal/microchunk_pass_1_3.simple'
  test_samples(10)%stkfile  = 'rejection_benchmarks/Betagal/microchunk_pass_1_3.mrc'
  test_samples(10)%msk      = 200

  test_samples(11)%projfile = 'rejection_benchmarks/Betagal/microchunk_pass_2_1.simple'
  test_samples(11)%stkfile  = 'rejection_benchmarks/Betagal/microchunk_pass_2_1.mrc'
  test_samples(11)%msk      = 200

  test_samples(12)%projfile = 'rejection_benchmarks/Betagal/microchunk_pass_2_2.simple'
  test_samples(12)%stkfile  = 'rejection_benchmarks/Betagal/microchunk_pass_2_2.mrc'
  test_samples(12)%msk      = 200

  test_samples(13)%projfile = 'rejection_benchmarks/Betagal/microchunk_pass_2_3.simple'
  test_samples(13)%stkfile  = 'rejection_benchmarks/Betagal/microchunk_pass_2_3.mrc'
  test_samples(13)%msk      = 200

  test_samples(14)%projfile = 'rejection_benchmarks/Betagal/microchunk_match_1.simple'
  test_samples(14)%stkfile  = 'rejection_benchmarks/Betagal/microchunk_match_1.mrc'
  test_samples(14)%msk      = 200

  test_samples(15)%projfile = 'rejection_benchmarks/Betagal/microchunk_match_2.simple'
  test_samples(15)%stkfile  = 'rejection_benchmarks/Betagal/microchunk_match_2.mrc'
  test_samples(15)%msk      = 200

  test_samples(16)%projfile = 'rejection_benchmarks/Betagal/microchunk_match_3.simple'
  test_samples(16)%stkfile  = 'rejection_benchmarks/Betagal/microchunk_match_3.mrc'
  test_samples(16)%msk      = 200

  test_samples(17)%projfile = 'rejection_benchmarks/Betagal/refchunk.simple'
  test_samples(17)%stkfile  = 'rejection_benchmarks/Betagal/refchunk.mrc'
  test_samples(17)%msk      = 200

  ! PepT2
  test_samples(18)%projfile = 'rejection_benchmarks/PepT2/microchunk_pass_1_1.simple'
  test_samples(18)%stkfile  = 'rejection_benchmarks/PepT2/microchunk_pass_1_1.mrc'
  test_samples(18)%msk      = 160

  test_samples(19)%projfile = 'rejection_benchmarks/PepT2/microchunk_pass_1_2.simple'
  test_samples(19)%stkfile  = 'rejection_benchmarks/PepT2/microchunk_pass_1_2.mrc'
  test_samples(19)%msk      = 160

  test_samples(20)%projfile = 'rejection_benchmarks/PepT2/microchunk_pass_1_3.simple'
  test_samples(20)%stkfile  = 'rejection_benchmarks/PepT2/microchunk_pass_1_3.mrc'
  test_samples(20)%msk      = 160

  test_samples(21)%projfile = 'rejection_benchmarks/PepT2/microchunk_pass_2_1.simple'
  test_samples(21)%stkfile  = 'rejection_benchmarks/PepT2/microchunk_pass_2_1.mrc'
  test_samples(21)%msk      = 160

  test_samples(22)%projfile = 'rejection_benchmarks/PepT2/microchunk_pass_2_2.simple'
  test_samples(22)%stkfile  = 'rejection_benchmarks/PepT2/microchunk_pass_2_2.mrc'
  test_samples(22)%msk      = 160

  test_samples(23)%projfile = 'rejection_benchmarks/PepT2/microchunk_pass_2_3.simple'
  test_samples(23)%stkfile  = 'rejection_benchmarks/PepT2/microchunk_pass_2_3.mrc'
  test_samples(23)%msk      = 160

  test_samples(24)%projfile = 'rejection_benchmarks/PepT2/microchunk_match_1.simple'
  test_samples(24)%stkfile  = 'rejection_benchmarks/PepT2/microchunk_match_1.mrc'
  test_samples(24)%msk      = 160

  test_samples(25)%projfile = 'rejection_benchmarks/PepT2/microchunk_match_2.simple'
  test_samples(25)%stkfile  = 'rejection_benchmarks/PepT2/microchunk_match_2.mrc'
  test_samples(25)%msk      = 160

  test_samples(26)%projfile = 'rejection_benchmarks/PepT2/microchunk_match_3.simple'
  test_samples(26)%stkfile  = 'rejection_benchmarks/PepT2/microchunk_match_3.mrc'
  test_samples(26)%msk      = 160

  test_samples(27)%projfile = 'rejection_benchmarks/PepT2/refchunk.simple'
  test_samples(27)%stkfile  = 'rejection_benchmarks/PepT2/refchunk.mrc'
  test_samples(27)%msk      = 160

  ! SLC
  test_samples(28)%projfile = 'rejection_benchmarks/SLC/microchunk_pass_1_1.simple'
  test_samples(28)%stkfile  = 'rejection_benchmarks/SLC/microchunk_pass_1_1.mrc'
  test_samples(28)%msk      = 150

  test_samples(29)%projfile = 'rejection_benchmarks/SLC/microchunk_pass_1_2.simple'
  test_samples(29)%stkfile  = 'rejection_benchmarks/SLC/microchunk_pass_1_2.mrc'
  test_samples(29)%msk      = 150

  test_samples(30)%projfile = 'rejection_benchmarks/SLC/microchunk_pass_1_3.simple'
  test_samples(30)%stkfile  = 'rejection_benchmarks/SLC/microchunk_pass_1_3.mrc'
  test_samples(30)%msk      = 150

  test_samples(31)%projfile = 'rejection_benchmarks/SLC/microchunk_pass_2_1.simple'
  test_samples(31)%stkfile  = 'rejection_benchmarks/SLC/microchunk_pass_2_1.mrc'
  test_samples(31)%msk      = 150

  test_samples(32)%projfile = 'rejection_benchmarks/SLC/microchunk_pass_2_2.simple'
  test_samples(32)%stkfile  = 'rejection_benchmarks/SLC/microchunk_pass_2_2.mrc'
  test_samples(32)%msk      = 150

  test_samples(33)%projfile = 'rejection_benchmarks/SLC/microchunk_pass_2_3.simple'
  test_samples(33)%stkfile  = 'rejection_benchmarks/SLC/microchunk_pass_2_3.mrc'
  test_samples(33)%msk      = 150

  test_samples(34)%projfile = 'rejection_benchmarks/SLC/microchunk_match_1.simple'
  test_samples(34)%stkfile  = 'rejection_benchmarks/SLC/microchunk_match_1.mrc'
  test_samples(34)%msk      = 150

  test_samples(35)%projfile = 'rejection_benchmarks/SLC/microchunk_match_2.simple'
  test_samples(35)%stkfile  = 'rejection_benchmarks/SLC/microchunk_match_2.mrc'
  test_samples(35)%msk      = 150

  test_samples(36)%projfile = 'rejection_benchmarks/SLC/microchunk_match_3.simple'
  test_samples(36)%stkfile  = 'rejection_benchmarks/SLC/microchunk_match_3.mrc'
  test_samples(36)%msk      = 150

  test_samples(37)%projfile = 'rejection_benchmarks/SLC/refchunk.simple'
  test_samples(37)%stkfile  = 'rejection_benchmarks/SLC/refchunk.mrc'
  test_samples(37)%msk      = 150

  ! Not1
  test_samples(38)%projfile = 'rejection_benchmarks/Not1/microchunk_pass_1_1.simple'
  test_samples(38)%stkfile  = 'rejection_benchmarks/Not1/microchunk_pass_1_1.mrc'
  test_samples(38)%msk      = 156

  test_samples(39)%projfile = 'rejection_benchmarks/Not1/microchunk_pass_1_2.simple'
  test_samples(39)%stkfile  = 'rejection_benchmarks/Not1/microchunk_pass_1_2.mrc'
  test_samples(39)%msk      = 156

  test_samples(40)%projfile = 'rejection_benchmarks/Not1/microchunk_pass_1_3.simple'
  test_samples(40)%stkfile  = 'rejection_benchmarks/Not1/microchunk_pass_1_3.mrc'
  test_samples(40)%msk      = 156

  test_samples(41)%projfile = 'rejection_benchmarks/Not1/microchunk_pass_2_1.simple'
  test_samples(41)%stkfile  = 'rejection_benchmarks/Not1/microchunk_pass_2_1.mrc'
  test_samples(41)%msk      = 156

  test_samples(42)%projfile = 'rejection_benchmarks/Not1/microchunk_pass_2_2.simple'
  test_samples(42)%stkfile  = 'rejection_benchmarks/Not1/microchunk_pass_2_2.mrc'
  test_samples(42)%msk      = 156

  test_samples(43)%projfile = 'rejection_benchmarks/Not1/microchunk_pass_2_3.simple'
  test_samples(43)%stkfile  = 'rejection_benchmarks/Not1/microchunk_pass_2_3.mrc'
  test_samples(43)%msk      = 156

  test_samples(44)%projfile = 'rejection_benchmarks/Not1/microchunk_match_1.simple'
  test_samples(44)%stkfile  = 'rejection_benchmarks/Not1/microchunk_match_1.mrc'
  test_samples(44)%msk      = 156

  test_samples(45)%projfile = 'rejection_benchmarks/Not1/microchunk_match_2.simple'
  test_samples(45)%stkfile  = 'rejection_benchmarks/Not1/microchunk_match_2.mrc'
  test_samples(45)%msk      = 156

  test_samples(46)%projfile = 'rejection_benchmarks/Not1/microchunk_match_3.simple'
  test_samples(46)%stkfile  = 'rejection_benchmarks/Not1/microchunk_match_3.mrc'
  test_samples(46)%msk      = 156

  test_samples(47)%projfile = 'rejection_benchmarks/Not1/refchunk.simple'
  test_samples(47)%stkfile  = 'rejection_benchmarks/Not1/refchunk.mrc'
  test_samples(47)%msk      = 156

  ! MotAB
  test_samples(48)%projfile = 'rejection_benchmarks/MotAB/microchunk_pass_1_1.simple'
  test_samples(48)%stkfile  = 'rejection_benchmarks/MotAB/microchunk_pass_1_1.mrc'
  test_samples(48)%msk      = 250

  test_samples(49)%projfile = 'rejection_benchmarks/MotAB/microchunk_pass_1_2.simple'
  test_samples(49)%stkfile  = 'rejection_benchmarks/MotAB/microchunk_pass_1_2.mrc'
  test_samples(49)%msk      = 250

  test_samples(50)%projfile = 'rejection_benchmarks/MotAB/microchunk_pass_1_3.simple'
  test_samples(50)%stkfile  = 'rejection_benchmarks/MotAB/microchunk_pass_1_3.mrc'
  test_samples(50)%msk      = 250

  test_samples(51)%projfile = 'rejection_benchmarks/MotAB/microchunk_pass_2_1.simple'
  test_samples(51)%stkfile  = 'rejection_benchmarks/MotAB/microchunk_pass_2_1.mrc'
  test_samples(51)%msk      = 250

  test_samples(52)%projfile = 'rejection_benchmarks/MotAB/microchunk_pass_2_2.simple'
  test_samples(52)%stkfile  = 'rejection_benchmarks/MotAB/microchunk_pass_2_2.mrc'
  test_samples(52)%msk      = 250

  test_samples(53)%projfile = 'rejection_benchmarks/MotAB/microchunk_pass_2_3.simple'
  test_samples(53)%stkfile  = 'rejection_benchmarks/MotAB/microchunk_pass_2_3.mrc'
  test_samples(53)%msk      = 250

  test_samples(54)%projfile = 'rejection_benchmarks/MotAB/microchunk_match_1.simple'
  test_samples(54)%stkfile  = 'rejection_benchmarks/MotAB/microchunk_match_1.mrc'
  test_samples(54)%msk      = 250

  test_samples(55)%projfile = 'rejection_benchmarks/MotAB/microchunk_match_2.simple'
  test_samples(55)%stkfile  = 'rejection_benchmarks/MotAB/microchunk_match_2.mrc'
  test_samples(55)%msk      = 250

  test_samples(56)%projfile = 'rejection_benchmarks/MotAB/microchunk_match_3.simple'
  test_samples(56)%stkfile  = 'rejection_benchmarks/MotAB/microchunk_match_3.mrc'
  test_samples(56)%msk      = 250

  test_samples(57)%projfile = 'rejection_benchmarks/MotAB/refchunk.simple'
  test_samples(57)%stkfile  = 'rejection_benchmarks/MotAB/refchunk.mrc'
  test_samples(57)%msk      = 250

  ! CLC
  test_samples(58)%projfile = 'rejection_benchmarks/CLC/microchunk_pass_1_1.simple'
  test_samples(58)%stkfile  = 'rejection_benchmarks/CLC/microchunk_pass_1_1.mrc'
  test_samples(58)%msk      = 209

  test_samples(59)%projfile = 'rejection_benchmarks/CLC/microchunk_pass_1_2.simple'
  test_samples(59)%stkfile  = 'rejection_benchmarks/CLC/microchunk_pass_1_2.mrc'
  test_samples(59)%msk      = 209

  test_samples(60)%projfile = 'rejection_benchmarks/CLC/microchunk_pass_1_3.simple'
  test_samples(60)%stkfile  = 'rejection_benchmarks/CLC/microchunk_pass_1_3.mrc'
  test_samples(60)%msk      = 209

  test_samples(61)%projfile = 'rejection_benchmarks/CLC/microchunk_pass_2_1.simple'
  test_samples(61)%stkfile  = 'rejection_benchmarks/CLC/microchunk_pass_2_1.mrc'
  test_samples(61)%msk      = 209

  test_samples(62)%projfile = 'rejection_benchmarks/CLC/microchunk_pass_2_2.simple'
  test_samples(62)%stkfile  = 'rejection_benchmarks/CLC/microchunk_pass_2_2.mrc'
  test_samples(62)%msk      = 209

  test_samples(63)%projfile = 'rejection_benchmarks/CLC/microchunk_pass_2_3.simple'
  test_samples(63)%stkfile  = 'rejection_benchmarks/CLC/microchunk_pass_2_3.mrc'
  test_samples(63)%msk      = 209

  test_samples(64)%projfile = 'rejection_benchmarks/CLC/microchunk_match_1.simple'
  test_samples(64)%stkfile  = 'rejection_benchmarks/CLC/microchunk_match_1.mrc'
  test_samples(64)%msk      = 209

  test_samples(65)%projfile = 'rejection_benchmarks/CLC/microchunk_match_2.simple'
  test_samples(65)%stkfile  = 'rejection_benchmarks/CLC/microchunk_match_2.mrc'
  test_samples(65)%msk      = 209

  test_samples(66)%projfile = 'rejection_benchmarks/CLC/microchunk_match_3.simple'
  test_samples(66)%stkfile  = 'rejection_benchmarks/CLC/microchunk_match_3.mrc'
  test_samples(66)%msk      = 209

  test_samples(67)%projfile = 'rejection_benchmarks/CLC/refchunk.simple'
  test_samples(67)%stkfile  = 'rejection_benchmarks/CLC/refchunk.mrc'
  test_samples(67)%msk      = 209

  ! Ribosome
  test_samples(68)%projfile = 'rejection_benchmarks/Ribosome/microchunk_pass_1_1.simple'
  test_samples(68)%stkfile  = 'rejection_benchmarks/Ribosome/microchunk_pass_1_1.mrc'
  test_samples(68)%msk      = 362

  test_samples(69)%projfile = 'rejection_benchmarks/Ribosome/microchunk_pass_1_2.simple'
  test_samples(69)%stkfile  = 'rejection_benchmarks/Ribosome/microchunk_pass_1_2.mrc'
  test_samples(69)%msk      = 362

  test_samples(70)%projfile = 'rejection_benchmarks/Ribosome/microchunk_pass_1_3.simple'
  test_samples(70)%stkfile  = 'rejection_benchmarks/Ribosome/microchunk_pass_1_3.mrc'
  test_samples(70)%msk      = 362

  test_samples(71)%projfile = 'rejection_benchmarks/Ribosome/microchunk_pass_2_1.simple'
  test_samples(71)%stkfile  = 'rejection_benchmarks/Ribosome/microchunk_pass_2_1.mrc'
  test_samples(71)%msk      = 362

  test_samples(72)%projfile = 'rejection_benchmarks/Ribosome/microchunk_pass_2_2.simple'
  test_samples(72)%stkfile  = 'rejection_benchmarks/Ribosome/microchunk_pass_2_2.mrc'
  test_samples(72)%msk      = 362

  test_samples(73)%projfile = 'rejection_benchmarks/Ribosome/microchunk_pass_2_3.simple'
  test_samples(73)%stkfile  = 'rejection_benchmarks/Ribosome/microchunk_pass_2_3.mrc'
  test_samples(73)%msk      = 362

  test_samples(74)%projfile = 'rejection_benchmarks/Ribosome/microchunk_match_1.simple'
  test_samples(74)%stkfile  = 'rejection_benchmarks/Ribosome/microchunk_match_1.mrc'
  test_samples(74)%msk      = 362

  test_samples(75)%projfile = 'rejection_benchmarks/Ribosome/microchunk_match_2.simple'
  test_samples(75)%stkfile  = 'rejection_benchmarks/Ribosome/microchunk_match_2.mrc'
  test_samples(75)%msk      = 362

  test_samples(76)%projfile = 'rejection_benchmarks/Ribosome/microchunk_match_3.simple'
  test_samples(76)%stkfile  = 'rejection_benchmarks/Ribosome/microchunk_match_3.mrc'
  test_samples(76)%msk      = 362

  test_samples(77)%projfile = 'rejection_benchmarks/Ribosome/refchunk.simple'
  test_samples(77)%stkfile  = 'rejection_benchmarks/Ribosome/refchunk.mrc'
  test_samples(77)%msk      = 362

  ! SusCDE
  test_samples(78)%projfile = 'rejection_benchmarks/SusCDE/microchunk_pass_1_1.simple'
  test_samples(78)%stkfile  = 'rejection_benchmarks/SusCDE/microchunk_pass_1_1.mrc'
  test_samples(78)%msk      = 190

  test_samples(79)%projfile = 'rejection_benchmarks/SusCDE/microchunk_pass_1_2.simple'
  test_samples(79)%stkfile  = 'rejection_benchmarks/SusCDE/microchunk_pass_1_2.mrc'
  test_samples(79)%msk      = 190

  test_samples(80)%projfile = 'rejection_benchmarks/SusCDE/microchunk_pass_1_3.simple'
  test_samples(80)%stkfile  = 'rejection_benchmarks/SusCDE/microchunk_pass_1_3.mrc'
  test_samples(80)%msk      = 190

  test_samples(81)%projfile = 'rejection_benchmarks/SusCDE/microchunk_pass_2_1.simple'
  test_samples(81)%stkfile  = 'rejection_benchmarks/SusCDE/microchunk_pass_2_1.mrc'
  test_samples(81)%msk      = 190

  test_samples(82)%projfile = 'rejection_benchmarks/SusCDE/microchunk_pass_2_2.simple'
  test_samples(82)%stkfile  = 'rejection_benchmarks/SusCDE/microchunk_pass_2_2.mrc'
  test_samples(82)%msk      = 190

  test_samples(83)%projfile = 'rejection_benchmarks/SusCDE/microchunk_pass_2_3.simple'
  test_samples(83)%stkfile  = 'rejection_benchmarks/SusCDE/microchunk_pass_2_3.mrc'
  test_samples(83)%msk      = 190

  test_samples(84)%projfile = 'rejection_benchmarks/SusCDE/microchunk_match_1.simple'
  test_samples(84)%stkfile  = 'rejection_benchmarks/SusCDE/microchunk_match_1.mrc'
  test_samples(84)%msk      = 190

  test_samples(85)%projfile = 'rejection_benchmarks/SusCDE/microchunk_match_2.simple'
  test_samples(85)%stkfile  = 'rejection_benchmarks/SusCDE/microchunk_match_2.mrc'
  test_samples(85)%msk      = 190

  test_samples(86)%projfile = 'rejection_benchmarks/SusCDE/microchunk_match_3.simple'
  test_samples(86)%stkfile  = 'rejection_benchmarks/SusCDE/microchunk_match_3.mrc'
  test_samples(86)%msk      = 190

  test_samples(87)%projfile = 'rejection_benchmarks/SusCDE/refchunk.simple'
  test_samples(87)%stkfile  = 'rejection_benchmarks/SusCDE/refchunk.mrc'
  test_samples(87)%msk      = 190

  do isample=1, size(test_samples)
    write(*,*) test_samples(isample)%projfile%to_char()
    call spproj%read(test_samples(isample)%projfile)
    call spproj%add_cavgs2os_out(test_samples(isample)%stkfile, spproj%get_smpd())
    cavg_imgs = read_cavgs_into_imgarr(spproj)

    call rejector%new(cavg_imgs, real(test_samples(isample)%msk))

    ! mirror rejection in microchunking
    if( test_samples(isample)%projfile%has_substr(string('refchunk')) ) then
      call rejector%reject_pop(spproj%os_cls2D, thres=0.0025)
    else if( test_samples(isample)%projfile%has_substr(string('microchunk_match')) ) then
      call rejector%reject_pop(spproj%os_cls2D, thres=0.0025)
    else if( test_samples(isample)%projfile%has_substr(string('microchunk_pass_2')) ) then
      call rejector%reject_pop(spproj%os_cls2D, thres=0.0035)
    else
      call rejector%reject_pop(spproj%os_cls2D)
    end if
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
    ! mirror rejection in microchunking
    if( test_samples(isample)%projfile%has_substr(string('refchunk')) ) then
      call rejector%reject_local_variance(strong_thresh=-2.0, weak_thresh=-2.0)
    else if( test_samples(isample)%projfile%has_substr(string('microchunk_match')) ) then
      call rejector%reject_local_variance(strong_thresh=-2.0, weak_thresh=-2.0)
    else if( test_samples(isample)%projfile%has_substr(string('microchunk_pass_2')) ) then
      call rejector%reject_local_variance(strong_thresh=-1.0, weak_thresh=-1.0)
    else
      call rejector%reject_local_variance()
    end if
    states = rejector%get_states()
    call spproj%os_cls2D%set_all('state', states)
    jpeg_out = swap_suffix(test_samples(isample)%stkfile, string('_local_variance_rejected.jpg'), string(MRC_EXT))
    call spproj%cavgs2jpg(cavg_inds, jpeg_out, xtiles, ytiles, ignore_states=.false., invert_states=.true.)

    cavg_imgs = read_cavgs_into_imgarr(spproj)
    call rejector%new(cavg_imgs, real(test_samples(isample)%msk))
    if( test_samples(isample)%projfile%has_substr(string('refchunk')) ) then
      call rejector%reject_pop(spproj%os_cls2D, thres=0.0025)
    else if( test_samples(isample)%projfile%has_substr(string('microchunk_match')) ) then
      call rejector%reject_pop(spproj%os_cls2D, thres=0.0025)
    else if( test_samples(isample)%projfile%has_substr(string('microchunk_pass_2')) ) then
      call rejector%reject_pop(spproj%os_cls2D, thres=0.0035)
    else
      call rejector%reject_pop(spproj%os_cls2D)
    end if
    call rejector%reject_res(spproj%os_cls2D)
    call rejector%reject_mask()
    if( test_samples(isample)%projfile%has_substr(string('refchunk')) ) then
      call rejector%reject_local_variance(strong_thresh=-2.0, weak_thresh=-2.0)
    else if( test_samples(isample)%projfile%has_substr(string('microchunk_match')) ) then
      call rejector%reject_local_variance(strong_thresh=-2.0, weak_thresh=-2.0)
    else if( test_samples(isample)%projfile%has_substr(string('microchunk_pass_2')) ) then
      call rejector%reject_local_variance(strong_thresh=-1.0, weak_thresh=-1.0)
    else
      call rejector%reject_local_variance()
    end if
    states = rejector%get_states()
    call spproj%os_cls2D%set_all('state', states)
    jpeg_out = swap_suffix(test_samples(isample)%stkfile, string('_selected.jpg'), string(MRC_EXT))
    call spproj%cavgs2jpg(cavg_inds, jpeg_out, xtiles, ytiles, ignore_states=.false., invert_states=.false.)

    call rejector%kill()
    call spproj%kill()
    call dealloc_imgarr(cavg_imgs)

  end do

end program main
