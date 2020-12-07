
module newton_mod
    use SLAU_solver
    use precision
    contains

    subroutine newton(X0,iter_count,X,f)
    real(mp), dimension(1:), intent(in) :: X0
    real(mp), dimension(1:size(X0)), intent(out) :: X
    real(mp), dimension(1:size(X0)) :: Xnew 
    real(mp), dimension(1:size(X),1:size(X)) :: df 
    real(mp), dimension(1:size(X)+1,1:size(X)) :: M 
    real(mp) :: eps
    integer, intent(in) :: iter_count
    integer :: i, n
    interface
        function f(X)
        use precision
        real(mp), dimension(1:), intent(in) :: X
        real(mp), dimension(1:size(X)) :: f
        end function f
    end interface
    n=size(X0)
    eps = 0.001_mp
    X=X0; i=iter_count; Xnew=X0+1.0_mp
    do while (sum(abs(X-Xnew))>eps .and. i>=1)
        X=Xnew
        n=size(X)
        call Jacobian(X,df,f)
        M(n+1,1:n)=-f(X)
        M(1:n,1:n)=df(1:n,1:n)
        call lead(M,Xnew)
        Xnew=Xnew+X
        i=i-1
    enddo
    end subroutine newton

    subroutine Jacobian(X0,df,f)
    implicit none
    real(mp), dimension(1:), intent(in) :: X0
    real(mp), dimension(1:size(X0),1:size(X0)), intent(out) :: df
    real(mp), dimension(1:size(X0)) :: X
    real(mp) :: eps!=0.01_mp**(mp+1)
    integer :: i, n
    interface
        function f(X)
        use precision
        real(mp), dimension(1:), intent(in) :: X
        real(mp), dimension(1:size(X)) :: f
        end function f
    end interface
    n=size(X0)
    eps = 0.001_mp
    do i=1,n
        X=X0
        X(i)=X0(i)+eps
        df(i,1:n)=(f(X)-f(X0))/eps
    enddo
    end subroutine Jacobian
   
 end module newton_mod