module potentials
    use precision
    use constants
    implicit none
    
    abstract interface
        function potential (x,y,z)
            real(8), intent(in) :: x,y,z
            real(8) :: potential
        end function potential
    end interface

    type function_pointer
        procedure (potential), pointer, nopass :: this_potential
    end type
    
    contains
    
    function Plummer_buldge_potential(x,y,z)
    real(mp), intent(in) :: x, y, z
    real(mp) :: Plummer_buldge_potential
    real(mp) :: R
    R = sqrt(x**2 + y**2)
    Plummer_buldge_potential = - M_buldge/sqrt(R**2 + b_buldge**2)
    end function Plummer_buldge_potential
    
    function Miyamoto_Nagai_buldge_potential(x,y,z)
    real(mp), intent(in) :: x, y, z
    real(mp) :: Miyamoto_Nagai_buldge_potential
    real(mp) :: R
    R = sqrt(x**2 + y**2)
    Miyamoto_Nagai_buldge_potential = M_buldge/sqrt(R**2 + (a_buldge+sqrt(z**2+b_buldge**2))**2)
    end function Miyamoto_Nagai_buldge_potential

    function Miyamoto_Nagai_disk_potential(x,y,z)
    real(mp), intent(in) :: x, y, z
    real(mp) :: Miyamoto_Nagai_disk_potential
    real(mp) :: R
    R = sqrt(x**2 + y**2)
    Miyamoto_Nagai_disk_potential = M_disk/sqrt(R**2 + (a_disk+sqrt(z**2+b_disk**2))**2)
    end function Miyamoto_Nagai_disk_potential
    
    function Miyamoto_Nagai_thin_disk_potential(x,y,z)
    real(mp), intent(in) :: x, y, z
    real(mp) :: Miyamoto_Nagai_thin_disk_potential
    real(mp) :: R
    R = sqrt(x**2 + y**2)
    Miyamoto_Nagai_thin_disk_potential = M_thin_disk/sqrt(R**2 + (a_thin_disk+sqrt(z**2+b_thin_disk**2))**2)
    end function Miyamoto_Nagai_thin_disk_potential
    
    function Miyamoto_Nagai_thick_disk_potential(x,y,z)
    real(mp), intent(in) :: x, y, z
    real(mp) :: Miyamoto_Nagai_thick_disk_potential
    real(mp) :: R
    R = sqrt(x**2 + y**2)
    Miyamoto_Nagai_thick_disk_potential = M_thick_disk/sqrt(R**2 + (a_thick_disk+sqrt(z**2+b_thick_disk**2))**2)
    end function Miyamoto_Nagai_thick_disk_potential

    function Logarithmic_halo_potential(x,y,z)
    real(mp), intent(in) :: x, y, z
    real(mp) :: Logarithmic_halo_potential
    real(mp) :: R, q
    R = sqrt(x**2 + y**2)
    q = 1.0 ! model parameter, it variates in range from 0 to 1
    Logarithmic_halo_potential = -v0**2*log(R**2 + a_halo**2 + z**2/q**2)
    end function Logarithmic_halo_potential
    
end module potentials