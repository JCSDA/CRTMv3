! File: CloudCoeff_Generate_default_format.F90
program CloudCoeff_Generate_default_format
  use Message_Handler
  use Type_Kinds,        only: Double
  use CloudCoeff_Define
  use CloudCoeff_IO
  use cloudcoeff_populate, only: fill_cloudcoeff_all
  implicit none

  integer, parameter :: wp = Double

  integer, parameter :: n_MW_Frequencies = 31
  integer, parameter :: n_IR_Frequencies = 61
  integer, parameter :: n_MW_Radii       = 10
  integer, parameter :: n_IR_Radii       = 10
  integer, parameter :: n_Temperatures   = 5
  integer, parameter :: n_Densities      = 3
  integer, parameter :: n_Legendre_Terms = 39
  integer, parameter :: n_Phase_Elements = 1

  type(CloudCoeff_type) :: CC
  integer :: s, i

  ! Create with the core dims only (I1,I2,I3,I4,I5,I6)
  call CloudCoeff_Create(CC, n_MW_Frequencies, n_MW_Radii, &
                             n_IR_Frequencies, n_IR_Radii, &
                             n_Temperatures,   n_Densities, n_Legendre_Terms, n_Phase_Elements)

  ! Set series/phase counts (and clamp if needed)
  CC%Max_Legendre_Terms = n_Legendre_Terms
  CC%n_Legendre_Terms   = n_Legendre_Terms
  CC%Max_Phase_Elements = n_Phase_Elements
  CC%n_Phase_Elements   = n_Phase_Elements

  ! Coordinate vectors (use _Double literals)
  CC%Frequency_MW = [ 10.65_Double, 18.7_Double, 23.8_Double, 31.4_Double, 36.5_Double, &
                       50.3_Double,  52.8_Double, 53.6_Double, 54.4_Double, 55.5_Double, &
                       57.29_Double, 88.2_Double, 89.0_Double, 91.7_Double,118.75_Double, &
                      150.0_Double, 157.0_Double,166.0_Double,167.0_Double,170.0_Double,  &
                      183.31_Double,190.0_Double,220.0_Double,314.0_Double,325.0_Double,  &
                      340.0_Double, 360.0_Double,380.0_Double,424.0_Double,448.0_Double,  &
                      664.0_Double ]

  do i = 1, n_IR_Frequencies
    CC%Frequency_IR(i) = 200.0_Double + (i-1)* (37037.0_Double - 200.0_Double) / real(n_IR_Frequencies-1, Double)
  end do

  CC%Reff_MW     = [ 5._Double, 15._Double, 30._Double, 50._Double,  80._Double, &
                     150._Double, 300._Double, 500._Double, 800._Double, 1500._Double ]

  CC%Reff_IR     = [ 2._Double, 3._Double, 4._Double, 5._Double, 6._Double, 8._Double, &
                    10._Double,15._Double,20._Double,50._Double ]

  CC%Temperature = [ 263.16_Double, 273.16_Double, 282._Double, 290._Double, 300._Double ]

  CC%Density     = [ 100._Double, 400._Double, 900._Double ]

  ! Compute and populate all arrays
  call fill_cloudcoeff_all(CC, 0)
  call CloudCoeff_Inspect(CC)
  stop 'bp'
  ! Write NetCDF (signature: filename first, then object)
  s = CloudCoeff_WriteFile('CloudCoeff.nc4', CC)
  if (s /= SUCCESS) stop 1

  call CloudCoeff_Destroy(CC)
end program CloudCoeff_Generate_default_format

