program simple_h5_ncem
  use simple_hdf5
  implicit none
  character(len=41), parameter :: filename    = "/processing/ercius/k2data/GLCAuJon30_NP1.hdf" ! file name
  character(len=3),  parameter :: groupmdf    = "mdf"    ! mdf group name
  character(len=6),  parameter :: groupimages = "images" ! images group name
  character(len=1),  parameter :: group0      = "0"      ! image set 0 group name
  character(len=5),  parameter :: dsetimage   = "image"  ! image set data set name (same for all data sets)
  integer(hid_t) :: file_id    ! file identifier
  integer(hid_t) :: mdf_id     ! mdf group identifier
  integer(hid_t) :: images_id  ! images group identifier
  integer(hid_t) :: g0_id      ! image set 0 group identifier
  integer(hid_t) :: image0_id  ! image set 0 data set identifier
  integer(hid_t) :: dataspace0 ! file dataspace identifier
  integer(hid_t) :: memspace   ! memory dataspace identifier
  integer, dimension(400,300,300) :: dset_full ! full data in memory
  integer, dimension(300,300) :: data1         ! hold only one image in memory
  integer(hsize_t), dimension(3) :: data_dimsfull = (/400,300,300/) ! full data set dimensions
  integer(hsize_t), dimension(2) :: data_dims1    = (/300,300/)     ! single image dimensions
  integer(hsize_t), dimension(3) :: offsetfull    = (/0,0,0/)
  integer(hsize_t), dimension(2) :: offset1       = (/0,0/)
  integer(hsize_t), dimension(3) :: countfull     = (/400,300,300/)
  integer(hsize_t), dimension(2) :: memcount1     = (/300,300/)
  integer(hsize_t), dimension(3) :: datacount1    = (/1,300,300/)
  integer :: imagenum, error

  ! initialize fortran interface.
  call h5open_f(error)
  
  ! open the file as read only
  call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, error) !h5f_acc_rdwr_f

  ! open the top mdf group
  call h5gopen_f(file_id,groupmdf,mdf_id,error)

  ! open the images group
  call h5gopen_f(mdf_id,groupimages,images_id,error)

  ! open the image set 0 group
  call h5gopen_f(images_id,group0,g0_id,error)

  ! open the data set
  call h5dopen_f(g0_id, dsetimage, image0_id, error)

  ! read the full 400x300x300 data set into memory
  call h5dread_f(image0_id, h5t_native_integer, dset_full, data_dimsfull, error)

  ! close the data set
  call h5dclose_f(image0_id, error)

  print *, 'first intensity value should be = to 524:'
  print *, (dset_full(1,1,1)) !the first value in glcaujon30_np1.hdf is 524

  ! now try to read only the first 1x300x300 part of the data
  
  ! open the data set
  call h5dopen_f(g0_id, dsetimage, image0_id, error)

  ! get the data space from the file  
  call h5dget_space_f(image0_id,dataspace0,error)
  
  ! select part of the data in the file dataset to read
  ! change offset to select different images
  imagenum = 5
  offsetfull = (/imagenum,0,0/)
  call h5sselect_hyperslab_f(dataspace0,h5s_select_set_f,offsetfull,datacount1,error)

  ! create space in memory for data to read into
  call h5screate_simple_f(2, data_dims1, memspace, error)

  ! select space in memory.
  call h5sselect_hyperslab_f(memspace, h5s_select_set_f, offset1, memcount1, error)

  call h5dread_f(image0_id, h5t_native_integer, data1, data_dims1, error, memspace, dataspace0)

  print *, 'these should match:'
  print *, '  single image read (1,1):'
  print *,data1(1,1)
  print *, '  full data read of same pixel:'
  print *, (dset_full(imagenum+1,1,1))

  ! close the data set
  call h5dclose_f(image0_id, error)

  ! close the groups
  call h5gclose_f(g0_id, error)
  call h5gclose_f(images_id, error)
  call h5gclose_f(mdf_id, error)

  ! terminate access to the file.
  call h5fclose_f(file_id, error)

  ! close fortran interface.
  call h5close_f(error)

end program simple_h5_ncem