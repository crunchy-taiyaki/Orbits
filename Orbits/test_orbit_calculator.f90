program cluster_orbit
    use precision
    use orbit_calculator
    implicit none
    real(mp) :: t_start, t_end
    integer :: t_point_number
    real(mp) :: init_conditions(0:1) ! x0, y0, z0, x0', y0', z0'
    character :: model
    real(mp), allocatable :: t(:), orbit_solution(:,:)
    integer :: i !counter
    
    model = '1'
    t_start = 0.0_mp
    t_end = 1.0_mp
    t_point_number = 10 !should be more 5
    init_conditions = (/-1.0, 0.0/)
    allocate(orbit_solution(0:t_point_number-1,0:size(init_conditions)-1))
    call calc_orbit(model, t_start, t_end, t_point_number, init_conditions, t, orbit_solution)
      
    do i=0, t_point_number-1
        write(*,*) t(i), orbit_solution(i,0) - (t(i)-1.0_mp)*exp(t(i))
    enddo
    deallocate(orbit_solution)
end program   