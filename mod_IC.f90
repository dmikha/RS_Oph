module ic
  use utils, only : dp
  implicit none
  private
!public objects
  !cross sections
  public interaction_angle_averaged_ic_cross_section_function
  public ic_losses_single_photon
contains


    elemental function interaction_angle_averaged_ic_cross_section_function(electron_energy,soft_photon_energy,&
         &gamma_energy) result(cross_section)
      ! C - all energies in mcc units
      use utils, only : pi, re=>electron_radius
      implicit none
      real(dp), intent(in) ::  electron_energy,gamma_energy,soft_photon_energy
      real(dp) cross_section
      real(dp) :: Emax,Emin,t1,t3,t7,t8,t9,t12,t17,t25,t34

      cross_section=0.
      Emax=electron_energy*(1._dp-0.25_dp/electron_energy**2)/(1._dp+0.25_dp/electron_energy/soft_photon_energy)!-1.d-2
      Emin=soft_photon_energy!+1.d-2
      cross_section=0._dp
      if( gamma_energy < Emax .and. gamma_energy > Emin) then
         t1 = electron_energy-gamma_energy
         t3 = 1._dp/gamma_energy
         t7 = log(4._dp*t1*soft_photon_energy*t3*electron_energy)
         t8 = electron_energy**2
         t9 = soft_photon_energy*t8
         t12 = soft_photon_energy*electron_energy*gamma_energy
         t17 = gamma_energy**2
         t25 = 1._dp/t1
         t34 = soft_photon_energy**2
         cross_section = pi*re**2*(-t7-(-4._dp*t9+4._dp*t12+gamma_energy)*(2._dp*t9-2._dp*t12+soft_photon_energy*t17+gamma_energy)*&
              &t3/electron_energy/soft_photon_energy*t25/4._dp)/t8/electron_energy*t25*gamma_energy/t34
      endif

      return
    end function interaction_angle_averaged_ic_cross_section_function
    


    elemental function ic_losses_single_photon(electron_energy,photon_energy) result(electron_losses)!per second!!!
      implicit none
      real(dp), intent(in) :: electron_energy,photon_energy
      real(dp) :: electron_losses

      real(dp) :: tmp1,tmp2,tmp3

      tmp1=4.e0_dp*electron_energy*photon_energy
      tmp2=log(1.e0_dp+0.16e0_dp*tmp1)
      tmp3=tmp1**2
      electron_losses= -4.16166e-14_dp*tmp2/(1.e0_dp+0.139e1_dp*tmp1)*&
           &(1._dp-0.46e-1_dp*tmp1/(1._dp+0.49e-1_dp*tmp3))*electron_energy
      return
    end function ic_losses_single_photon
    

end module ic
