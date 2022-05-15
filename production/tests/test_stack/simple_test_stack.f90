program simple_test_stack
    include 'simple_lib.f08'
    use simple_stack
    implicit none
    type(stack) :: stk
    integer     :: val
    call stk%new
    call stk%push(5)
    call stk%push(6)
    call stk%push(10)
    write(*, *) stk%contains(1), stk%contains(5), stk%contains(10)
    val = stk%pop()
    write(*, *) stk%contains(1), stk%contains(5), stk%contains(10)
    val = stk%pop()
    val = stk%pop()
    write(*, *) stk%contains(1), stk%contains(5), stk%contains(10)
    val = stk%pop()
    write(*, *) val, stk%is_empty()
end program simple_test_stack