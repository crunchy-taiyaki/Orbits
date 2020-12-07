module orbit_calculator
    use precision
    use constants
    use potentials
    use ODE_solver
    implicit none
    contains
    
function Phi_model_1(x,y,z) result(Phi_res)!Galactic potential: Plummer buldge, MN disk and logarithmic halo
    real(mp), intent(in) :: x, y, z
    real(mp) :: R
    real(mp) :: Phi_res
    R = sqrt(x**2 + y**2)
    Phi_res = Plummer_buldge_potential(x,y,z)+Miyamoto_Nagai_disk_potential(x,y,z)+Logarithmic_halo_potential(x,y,z)
end function Phi_model_1

function Phi_model_2(x,y,z) result(Phi_res)!Galactic potential: MN buldge, MN disk and logarithmic halo
    real(mp), intent(in) :: x, y, z
    real(mp) :: R
    real(mp) :: Phi_res
    R = sqrt(x**2 + y**2)
    Phi_res = Miyamoto_Nagai_buldge_potential(x,y,z)+Miyamoto_Nagai_disk_potential(x,y,z)+Logarithmic_halo_potential(x,y,z)
end function Phi_model_2

function Phi_model_3(x,y,z) result(Phi_res)!Galactic potential: Plummer buldge, MN thin&thick disks and logarithmic halo
    real(mp), intent(in) :: x, y, z
    real(mp) :: R
    real(mp) :: Phi_res
    R = sqrt(x**2 + y**2)
    Phi_res = Plummer_buldge_potential(x,y,z)+&
          &Miyamoto_Nagai_thin_disk_potential(x,y,z)+Miyamoto_Nagai_thick_disk_potential(x,y,z)+&
          &Logarithmic_halo_potential(x,y,z)
end function Phi_model_3

function differential_eq(t,X)
    real(mp), intent(in) :: t, X(0:)
    real(mp) :: differential_eq(0:size(X)-1)
    real(mp) :: eps
    
    !common block
    type(function_pointer) :: Phi
    common /galactic_potential/ Phi
    eps = 1.0d-3
    differential_eq(0) = X(3)
    differential_eq(1) = X(4)
    differential_eq(2) = X(5)
    differential_eq(3) = -Phi%this_potential(X(0)+eps, X(1), X(2))-Phi%this_potential(X(0)-eps, X(1), X(2))/(2.0_mp*eps)
    differential_eq(4) = -Phi%this_potential(X(0), X(1)+eps, X(2))-Phi%this_potential(X(0), X(1)-eps, X(2))/(2.0_mp*eps)
    differential_eq(5) = -Phi%this_potential(X(0), X(1), X(2)+eps)-Phi%this_potential(X(0), X(1), X(2)-eps)/(2.0_mp*eps)    
end function differential_eq

function test_fun(t,X)
    real(mp), intent(in) :: t, X(0:)
    real(mp) :: test_fun(0:size(X)-1)
    test_fun(0) = X(0) + exp(t)
    test_fun(1) = X(1) + exp(t)
end function test_fun

subroutine time_grid(t_start, t_end, t_point_number,t,t_step)
    real(mp), intent(in) :: t_start, t_end
    integer, intent(in) :: t_point_number
    real(mp), intent(inout) :: t(0:), t_step
    integer :: i !counter
    t_step = (t_end-t_start)/(t_point_number-1)
    do i=0,t_point_number-1
        t(i) = t_start + i*t_step
    enddo
end subroutine

subroutine calc_orbit(model, t_start, t_end, t_point_number, init_conditions, t, solution)
    character, intent(in) :: model
    real(mp), intent(in) :: t_start, t_end
    integer, intent(in) :: t_point_number
    real(mp), intent(in) :: init_conditions(0:)
    real(mp), intent(out) :: solution(0:t_point_number-1,0:size(init_conditions)-1)
    real(mp), allocatable, intent(out) :: t(:)
    real(mp) :: t_step
    integer :: i! counter
    
    !common block
    type(function_pointer) :: Phi
    common /galactic_potential/ Phi
    
    allocate(t(0:t_point_number))
    call time_grid(t_start, t_end, t_point_number, t, t_step)
    
    if (model.eq.'1') then
        Phi%this_potential => Phi_model_1
    else if (model.eq.'2') then
        Phi%this_potential => Phi_model_2
    else if (model.eq.'3') then
        Phi%this_potential => Phi_model_3
    else
        write(*,*)"Please, choose one from those: 1, 2 or 3"
    endif   
    !call interpolated_adams_solver(differential_eq,t,t_step,init_conditions,4,solution) !galactic potential should be defined in common block /galactic potential/
    call interpolated_adams_solver(test_fun,t,t_step,init_conditions,4,solution)
end subroutine calc_orbit

end module orbit_calculator