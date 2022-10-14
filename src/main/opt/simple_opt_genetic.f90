! simple probabilistic genetic algorithm
module simple_opt_genetic
    include 'simple_lib.f08'
    implicit none

    abstract interface
        function objective_func(bitstring) result(val)
            integer, intent(in)  :: bitstring(:)
            real                 :: val
        end function objective_func
    end interface
contains
    subroutine crossover(bitstr1, bitstr2, cross_rate)
        integer, intent(inout) :: bitstr1(:), bitstr2(:)
        real   , intent(in)    :: cross_rate
        integer :: cross_pt, N
        N = size(bitstr1)
        if( rand() < cross_rate )then
            cross_pt = 1 + floor(rand()*(N-1))
            bitstr2(         1:cross_pt) = bitstr1(         1:cross_pt)
            bitstr1(cross_pt+1:N)        = bitstr2(cross_pt+1:N)
        endif
    end subroutine crossover

    subroutine mutation(bitstring, nrefs, mut_rate, bits)
        use simple_stack
        integer,     intent(inout) :: bitstring(:)
        integer,     intent(in)    :: nrefs
        real   ,     intent(in)    :: mut_rate
        type(stack), intent(inout) :: bits
        integer     :: k
        real        :: val
        call bits%clear()
        do k = 1, nrefs
            call bits%push(k)
        enddo
        call bits%remove(real(bitstring(k)))
        do k = 1, size(bitstring)
            if( rand() < mut_rate )then
                bitstring(k) = bits%get_at(floor(rand()*(nrefs-1)))
            endif
        enddo
    end subroutine mutation

    function select_pop(pop, costs, times) result(ind)
        integer, intent(in)  :: pop(:,:)
        real,    intent(in)  :: costs(:)
        integer, intent(in)  :: times
        integer :: ind
        integer :: k, tmp_ind
        ind = 1+floor(rand()*size(pop,1))
        do k = 1, times
            tmp_ind = 1+floor(rand()*size(pop,1))
            if( costs(tmp_ind) > costs(ind) ) ind = tmp_ind
        enddo
    end function select_pop

    subroutine genetic_opt(obj_func, num_bits, nrefs, max_iter, pop_size, cross_rate, mut_rate, best, best_val, tol)
        use simple_stack
        procedure(objective_func), pointer :: obj_func
        integer, intent(in)    :: num_bits
        integer, intent(in)    :: nrefs
        integer, intent(in)    :: max_iter
        integer, intent(in)    :: pop_size
        real   , intent(in)    :: cross_rate
        real   , intent(in)    :: mut_rate
        integer, intent(inout) :: best(:)
        real   , intent(inout) :: best_val
        real   , intent(in), optional :: tol
        integer, allocatable   :: population(:,:), selected(:,:)
        real   , allocatable   :: costs(:)
        type(stack)            :: bits
        integer :: k,l
        real    :: prev_best_val
        call srand(time()+13)
        allocate(population(pop_size, num_bits), selected(pop_size, num_bits), source = 0)
        allocate(costs(pop_size), source = 0.)
        call bits%new(nrefs)
        do k = 1, pop_size
            do l = 1, num_bits
                population(k,l) = floor(rand()*nrefs)
            enddo
        enddo
        best          = population(1,:)
        best_val      = obj_func(best)
        prev_best_val = 0.
        do k = 1, max_iter
            do l = 1, pop_size
                costs(l) = obj_func(population(l,:))
                if( costs(l) > best_val )then
                    best          = population(l,:)
                    prev_best_val = best_val
                    best_val      = costs(l)
                endif
            enddo
            print *, 'current cost = ', best_val
            do l = 1, pop_size, 2
                selected(l  ,:) = population(select_pop(population, costs, times = 3), :)
                selected(l+1,:) = population(select_pop(population, costs, times = 3), :)
                call crossover(selected(l,:), selected(l+1,:), cross_rate)
                call mutation(selected(l  ,:), nrefs, mut_rate, bits)
                call mutation(selected(l+1,:), nrefs, mut_rate, bits)
            enddo
            population = selected
            if( best_val == num_bits ) return
            if( present(tol) .and. abs(prev_best_val - best_val) < tol ) return
        enddo
    end subroutine genetic_opt
end module simple_opt_genetic
