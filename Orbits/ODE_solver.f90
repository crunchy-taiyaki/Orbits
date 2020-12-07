module ODE_solver
    use precision
    use newton_mod
    contains    
    
    subroutine rungekut(this_function,X0,t,t_step,X)
    implicit none
    interface
        function this_function(t,x)
        real(8), intent(in) :: t,x(0:)
        real(8) :: this_function(0:size(x)-1)
        end function this_function
    end interface
    real(mp), dimension(1:), intent(in) :: X0
    real(mp), intent(in) :: t, t_step
    real(mp), dimension(1:size(X0)), intent(out) :: X
    real(mp), dimension(1:size(X0)) :: K1, K2, K3, K4

    K1=t_step*this_function(t,X0)
    K2=t_step*this_function(t+t_step/2,X0+K1/2)
    K3=t_step*this_function(t+t_step/2,X0+K2/2)
    K4=t_step*this_function(t+t_step,X0+K3)
    
    X=X0+(K1+2.0_mp*K2+2.0_mp*K3+K4)/6.0_mp
    end subroutine rungekut
    
        
    subroutine extra_adams(this_function,X0,t0,t_step,Xnew,X,adams_odder)    
    implicit none
    interface
        function this_function(t,x)
        real(8), intent(in) :: t,x(0:)
        real(8) :: this_function(0:size(x)-1)
        end function this_function
    end interface
    integer, intent(in) :: adams_odder
    real(mp), dimension(1:), intent(in) :: X0
    real(mp), intent(in) :: t0, t_step
    real(mp), dimension(1:size(X0)), intent(out) :: Xnew
    real(mp), dimension(0:adams_odder,1:size(X0)) :: X
    integer :: j, n
    n = adams_odder
    Xnew=X0
    do j=0,n-1
        Xnew=Xnew+A(n,j)*this_function(t0-t_step*(j-1),X(n-1-j,:))*t_step
    enddo
    X(n,:)=Xnew
    X(0:n-1,:)=X(1:n,:)
    
    end subroutine extra_adams
    
    subroutine interpolated_adams_solver(this_function,t,t_step,init_conditions,adams_odder,result_x)
    interface
        function this_function(t,x)
        real(8), intent(in) :: t,x(0:)
        real(8) :: this_function(0:size(x)-1)
        end function this_function
    end interface
    real(mp), intent(in) :: t(0:), t_step, init_conditions(0:)
    integer :: adams_odder
    real(mp), intent(out) :: result_x(0:size(t)-1,0:size(init_conditions)-1)
    real(mp), allocatable :: dummy(:), adams_interpolation_matrix(:,:) 
    integer :: i, j, t_size
    
    allocate(dummy(0:size(init_conditions)-1))
    allocate(adams_interpolation_matrix(-1:adams_odder-1,0:size(init_conditions)-1))

    adams_interpolation_matrix(-1,:) = init_conditions !t=0            
    do j = -1,(adams_odder-2)-1             
        call rungekut(this_function,adams_interpolation_matrix(j,:), t(j+1), t_step, adams_interpolation_matrix(j+1,:))
    enddo
    result_x(0:adams_odder-1,:) = adams_interpolation_matrix(-1:adams_odder-2,:)
    
    t_size = size(t)
    do i=adams_odder, t_size-1
        dummy = 0
        do j=0,adams_odder-2
            dummy = dummy + B(adams_odder,j)*this_function(t(i-j),adams_interpolation_matrix(i-j,:))
        enddo
        call newton(result_x(i-1,:),10,result_x(i,:),f_addams)
        adams_interpolation_matrix(-1:addams_odder-2,:)=result_x(i-(addams_odder-1):i,:)
    enddo            
    deallocate(adams_interpolation_matrix)
    
    contains
    function f_addams(x)
    implicit none
    real(mp), dimension(1:), intent(in) :: x(0:)
    real(mp), dimension(1:size(X)) :: f_addams
    f_addams = (B(adams_odder,-1)*this_function(t(i),x)+dummy)*t_step + result_x(i-1,:) - x 
    end function f_addams
    
    end subroutine interpolated_adams_solver
    
    elemental function A(n,j)
    implicit none
    integer, intent(in) :: n, j
    integer :: i
    real(mp) :: A
    real(mp) :: integral
    call integral_from_A(n,j,integral)
    A=(-1.0_mp)**j/product((/(i,i=1,j)/))/product((/(i,i=1,n-1-j)/))*integral
    
    contains
    
    pure subroutine integral_from_A(n,j,res)
    ! *** Процедура считает интеграл от функции z(z+1)..(z+n-1)/(z+j) в пределах [0,1]
    implicit none
    integer, intent(in) :: n, j
    integer :: a,i
    real(mp), dimension(0:n-1) :: z, helper
    real(mp),intent(out) :: res
    z=0;z(0)=1.0_mp;
    do a=0,n-1
        if(a/=j) then
            helper = z*a 
            z = eoshift(z, shift = -1) 
            z = z + helper
        endif
    enddo
    forall(i=0:n-1) z(i)=z(i)/(i+1)
    res = sum(z)    
    end subroutine integral_from_A
    
    end function A
    
    function B(n,j)
    implicit none
    integer, intent(in) :: n, j
    integer :: i
    real(mp) :: B
    real(mp) :: integral
    call integral_from_B(n,j,integral)
    if (j==-1) then
        B=1.0_mp/product((/(i,i=1,n-1)/))*integral
    else
        B=(-1.0_mp)**(j+1)/product((/(i,i=1,j+1)/))/product((/(i,i=1,n-2-j)/))*integral
    endif
    contains
    
    pure subroutine integral_from_B(n,j,res)
    ! *** Процедура считает интеграл от функции (z-1)z(z+1)..(z+n-2)/(z+j) в пределах [0,1]
    implicit none
    integer, intent(in) :: n, j
    integer :: a,i
    real(mp), dimension(-1:n-2) :: z, helper
    real(mp),intent(out) :: res
    z=0;z(-1)=1.0_mp;
    do a=-1,n-2  
        if(a/=j) then
            helper = z*a 
            z = eoshift(z, shift = -1) 
            z = z + helper
        endif
    enddo
    forall(i=-1:n-2) z(i)=z(i)/(i+2)
    res=sum(z)    
    end subroutine integral_from_B
    
    end function B
    
end module ODE_solver