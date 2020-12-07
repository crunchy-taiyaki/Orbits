module constants
    implicit none
    !buldge
    real(8), parameter :: pi=3.141592654
    real(8), parameter :: G = 4.300d-6 !gravitational constant in [kpc*Msun^-1*(km/s)^2]
    real(8), parameter :: v0 = 225.0 ! velocity in [km/s^-1]
    real(8), parameter :: M_buldge = 7.0 !buldge mass in [G*Msun]
    real(8), parameter :: a_buldge = 0.04 !buldge characteristic scale in [kpc]
    real(8), parameter :: b_buldge = 0.2 !buldge characteristic scale in [kpc]
    
    !disk
    real(8), parameter :: M_disk = 95 !disk mass in [G*Msun]
    real(8), parameter :: a_disk = 4.26 !disk characteristic scale in [kpc]
    real(8), parameter :: b_disk = 0.47 !disk characteristic scale in [kpc]
    
    !thin disk
    real(8), parameter :: M_thin_disk = 3.5d10*G !thin_disk mass in [G*Msun]
    real(8), parameter :: a_thin_disk = 2.6 !thin_disk characteristic scale in [kpc]
    real(8), parameter :: b_thin_disk = 0.3 !thin_disk characteristic scale in [kpc]
    
    !thick disk
    real(8), parameter :: M_thick_disk = 6.0d9*G !thick_disk mass in [G*Msun]
    real(8), parameter :: a_thick_disk = 2.0 !thick_disk characteristic scale in [kpc]
    real(8), parameter :: b_thick_disk = 0.9 !thick_disk characteristic scale in [kpc]
    
    !halo
    real(8), parameter :: a_halo =  3.20 !halo characteristic scale in [kpc]
end module constants