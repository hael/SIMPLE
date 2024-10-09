program simple_test_class_sample
include 'simple_lib.f08'
implicit none
type(class_sample), allocatable :: cs(:)
type(class_sample) :: cs_entry 
real, allocatable  :: rarr(:)
integer :: recsz, funit, io_stat
allocate(cs(3))
cs(1)%clsind = 1
cs(2)%clsind = 2
cs(3)%clsind = 3
cs(1)%pop    = 1
cs(2)%pop    = 2
cs(3)%pop    = 3
allocate(cs(1)%pinds(1), source=[1])
allocate(cs(2)%pinds(2), source=[1,2])
allocate(cs(3)%pinds(3), source=[1,2,3])
allocate(cs(1)%ccs(1),   source=[1.])
allocate(cs(2)%ccs(2),   source=[1.,0.])
allocate(cs(3)%ccs(3),   source=[1.,0.,-1.])

call print_class_sample(cs(1))
call print_class_sample(cs(2))
call print_class_sample(cs(3))

print *, '*********************'

call write_class_samples(cs, 'clssmp.bin')
call read_class_samples(cs, 'clssmp.bin')

call print_class_sample(cs(1))
call print_class_sample(cs(2))
call print_class_sample(cs(3))

end program simple_test_class_sample
