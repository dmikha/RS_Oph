module pair
  use utils, only : dp
  implicit none
  private
  public optical_depth_point_like_photon_source_array
contains

  function optical_depth_point_like_photon_source_array&
       &(gamma_energy,photon_energy,photon_density,initial_angle) result(depth)
    use utils, only : spacing=>linspace_middle_point,&
         &integral=>simple_integration
    implicit none
    real(dp), intent(in) :: gamma_energy(:),photon_energy(:),photon_density(:),initial_angle
    real(dp) :: depth(size(gamma_energy))

    real(dp) :: current_distance
    real(dp), parameter :: max_distance=20.
    integer, parameter :: N_steps=300

    real(dp) :: dtau(N_steps,size(gamma_energy)),distance(N_steps),cos_angles(N_steps),r_to_source(N_steps)
    real(dp) :: cross_section(size(photon_energy)),photon_density_at(size(photon_energy))

    integer :: i,j!,N_ph

    !N_ph=size(photon_energy)
    call spacing(distance,0._dp,max_distance)
    
    r_to_source=sqrt(1._dp+distance**2+2._dp*distance*cos(initial_angle))
    
    do i=1,N_steps
       if ( r_to_source(i) <= 1.e-3_dp ) stop "optical_depth_point_like_photon_source: r=0"
    enddo
    
    cos_angles=(distance+cos(initial_angle))/r_to_source

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) private(j,e_min,photon_energy,photon_density,cross_section)
    do i=1,size(gamma_energy)
       do j=1,N_steps
          if ( gamma_energy(i)*(1.-cos_angles(j)) > 0._dp) then
             !e_min=2._dp/(gamma_energy(i)*(1.-cos_angles(j)))
             !call logspacing(photon_energy,e_min,1.e5_dp*e_min)
             photon_density_at = photon_density / r_to_source(j)**2
             cross_section=total_pair_cross_section_function&
                  &(gamma_energy(i)*photon_energy*(1.-cos_angles(j)))
             dtau(j,i)=(1.-cos_angles(j))*integral(photon_energy,photon_density_at*cross_section)
          endif
       end do
       depth(i)=integral(distance,dtau(:,i))
    end do
!$OMP END PARALLEL DO
    return
  end function optical_depth_point_like_photon_source_array

  elemental function total_pair_cross_section_function&
       &(four_momentum_product) result(cross_section)
    use utils, only : pi, re=> electron_radius
    implicit none
    real(dp), intent(in) :: four_momentum_product
    real(dp) :: cross_section
    real(dp) :: t1,t2,t3
    cross_section=0._dp
    if ( four_momentum_product > 2._dp ) then
       t2=2._dp/four_momentum_product
       t1=sqrt(1._dp-t2)
       t3=(3._dp-t1**4)*(log(1._dp+t1)-log(1._dp-t1))
       cross_section= pi*re**2/2.*t2*(t3-2._dp*t1*(2._dp-t1**2))
       !pi*re**2/2.=1.25e-25_dp
    endif
    return
  end function total_pair_cross_section_function

  
end module pair


