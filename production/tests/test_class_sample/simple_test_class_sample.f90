program simple_test_class_sample
include 'simple_lib.f08'
implicit none
type(class_sample), allocatable :: cs(:), cs2(:)

call make_cs(cs)
call make_cs(cs2)

call print_class_sample(cs(1))
call print_class_sample(cs(2))
call print_class_sample(cs(3))

print *, '*********************'

call write_class_samples(cs, string('clssmp.bin'))
call read_class_samples(cs2, string('clssmp.bin'))

print *, class_samples_same(cs(1), cs2(1))
print *, class_samples_same(cs(2), cs2(2))
print *, class_samples_same(cs(3), cs2(3))

contains

    subroutine make_cs( cs )
        type(class_sample), allocatable, intent(inout) :: cs(:)
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
    end subroutine make_cs

end program simple_test_class_sample
