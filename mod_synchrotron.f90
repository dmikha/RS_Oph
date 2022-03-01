module synchrotron
  use utils, only : dp
  implicit none
  real(dp), parameter, public :: specific_critical_freq=3.398270086611665e-14_dp ! critical_freq/B-field strength {3e\hbar\over 2m_e^2c^3}, mcc/G units
  real(dp), parameter, public :: specific_energy_loss_rate=1.938485349573455e-09_dp! {2e**4\over 3m_e^3c^5} mcc/G**2/s units !
  real(dp), parameter, public :: specific_intensity_coef=35380.821101595771 ! mcc/G/s {\sqrt(3)\over 2\pi}{e^3\over mc^2\hbar} previously used 7.e4_dp -- 1.86 was accounted in the formulae
  private
!public objects
  public synchrotron_losses_function_iso !in isotropic field
  public synchrotron_emission_iso
contains

  pure function g(x,a) result(g_value)
    real(dp), intent(in) :: x,a(4)
    real(dp) :: g_value
    real(dp) :: tmp
    tmp=1._dp+a(3)*x**a(4)
    if ( tmp /= 0. ) then
       tmp=a(1)*x**a(2)/tmp+1._dp
       if (tmp /= 0. ) then 
          g_value=1./tmp
       else
          g_value=0.
       endif
    else
       g_value=0.
    endif
    return
  end function g

  
  elemental function synchrotron_losses_function_iso(electron_energy,B_field) result(electron_losses)
    implicit none
    real(dp), intent(in) :: electron_energy,B_field
    real(dp) :: electron_losses

    electron_losses = -specific_energy_loss_rate*(B_field*electron_energy)**2*0.66666!0.6666 comes from angle averaging
  end function synchrotron_losses_function_iso

  elemental function synchrotron_emission_kernel_iso(electron_energy,photon_energy,B_field) result(intensity) ! for chaotic B-fields
    implicit none
    real(dp), intent(in) :: electron_energy,photon_energy,B_field
    real(dp) :: intensity
    real(dp) :: critical_freq,intensity_coef
    real(dp) :: x,xx,y

    real(dp), dimension(4), parameter :: a_syn=[-0.552, 0.246, 0.856, 0.535] ! from Documents/iso_synchrotron.py

    critical_freq=specific_critical_freq*B_field*electron_energy**2
    intensity_coef=specific_intensity_coef*B_field
    x=photon_energy/critical_freq
    xx=x**0.33333
    y=1.808*xx/(1+1.151*xx)*exp(-x)*g(x,a_syn) ! from Documents/iso_synchrotron.py
    intensity=intensity_coef*y
    return
  end function synchrotron_emission_kernel_iso


  function synchrotron_emission_function_iso(electron_energy,electron_density,photon_energy,B_field) result(spectrum)
    use utils, only: integral=>simple_integration
    implicit none
    real(dp), intent(in) ::  electron_energy(:),electron_density(:),photon_energy(:),B_field
    real(dp) ::spectrum(size(photon_energy))
    
    integer nph
    real(dp) tmp

    if ( size(electron_energy) /= size(electron_density) ) &
         &stop "electron's arrays are inconcistent synchrotron_emission_function"

    do nph=1,size(photon_energy)
       tmp=integral(electron_energy,&
            &synchrotron_emission_kernel_iso(electron_energy,photon_energy(nph),B_field)*electron_density)
       spectrum(nph)=tmp/photon_energy(nph)
    enddo!
    return
  end function synchrotron_emission_function_iso

  subroutine synchrotron_emission_iso(electron_energy,electron_density,photon_energy,B_field,spectrum)
    implicit none
    real(dp), intent(in) ::  electron_energy(:),electron_density(:),photon_energy(:),B_field
    real(dp), intent(out)::  spectrum(:)
    if (size(spectrum) /= size(photon_energy) ) stop "array mismatch synchrotron_emission"
    spectrum=synchrotron_emission_function_iso(electron_energy,electron_density,photon_energy,B_field)
    return
  end subroutine synchrotron_emission_iso
  
end module synchrotron
  
