module units
use kind_parameters, only: sp, dp, i8

contains

!Energy conversion from keV to Eh (au)
real(dp) function EinAU(EinkeV)
implicit none
real(dp), intent(in) :: EinkeV
	EinAU = EinkeV*1.d3/27.21139
end function EinAU

!Distance conversion from angstroms to Bohr radii a0 (au)
real(dp) function dinAU(dinangstroms)
implicit none
real(dp), intent(in) :: dinangstroms
	dinAU = dinangstroms/0.5291772
end function dinAU

!Distance conversion from Bohr radii a0 (au) to meters (SI)
real(dp) function dinSI(dinAU)
implicit none
real(dp), intent(in) :: dinAU
	dinSI = dinAU*0.5291772d-10
end function dinSI

!Distance conversion from meters (SI) to Bohr radii a0 (au)
real(dp) function dinAUfromSI(dinSI)
implicit none
real(dp), intent(in) :: dinSI
	dinAUfromSI = dinSI/0.5291772d-10
end function dinAUfromSI

end module