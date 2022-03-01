module setup
  use utils, only : fp=>dp, AU=>AU_to_cm,pi,TeV_to_mcc
  use myparameters
  implicit none
  
  public physical_distance
  public photon_distance
  public shock_radius
  public set_up_photons
  public get_photons
  public get_Bfield_time
  public get_Bfield
  public get_density_time
  public get_density
  public beta_speed
contains
  function beta_speed(rrg) result(v)
    ! stellar wind speed
    ! v_w-> v_\infty
    implicit none
    real(fp), intent(in) :: rrg
    real(fp)             :: v
    real(fp), parameter  :: beta = 3
    v = v_w * ( 1 - R_s / rrg)**beta
  end function beta_speed
  
  function photon_distance(t) result(r)
    !distance from the explostion origin
    implicit none
    integer, intent(in) :: t
    real(fp)             :: r
    r = shock_radius(t) 
  end function photon_distance

  function shock_radius(t, vsh,vej) result(rsh)
    !shock radius and speed (and ejecta speed) calculation
    implicit none
    integer, intent(in)    :: t
    real(fp), intent(out), optional :: vsh,vej
    real(fp)                        :: rsh,v_ej
    real(fp)                        :: d_log_rho, v, rho1,rho2,rej
    integer                         :: l
    
    ! initial conditions: speed for a sperical explosion in constant density enviroment
    v    = v_ej0 * 1.3333 ** 0.333333
    v_ej = v_ej0
    rsh = (Time(2) - Time(1)) * v
    rej = (Time(2) - Time(1)) * v_ej 
    do l=2,min(t,N_t)
       rho1 = get_density(physical_distance(rsh     ))
       rho2 = get_density(physical_distance(rsh*1.01))
       !density log derivative, numerical estimate
       d_log_rho = (rho2 - rho1 ) / (rho1 * rsh *1.e-2)
       ! shock speed
       v = 4 * fill * rej**2 * v_ej /&
            & ( (4*fill-1) * rsh**2 + (rsh**3-rej**3) * d_log_rho * 1.33 * fill)
       ! ejecta slow dowing by the pressure
       v_ej = v_ej - 9.42 * rej**2 * rho1 * v**2 / m_ej
       ! new position of the ejecta and the shock
       rej = rej + (Time(l+1) - Time(l)) * v_ej
       rsh = rsh + (Time(l+1) - Time(l)) * v
    enddo
    if (present(vsh) ) then
       vsh = v
    endif
    if (present(vej) ) then
       vej = v_ej
    endif
  end function shock_radius
  
  elemental function physical_distance(r) result(r_cond)
    ! distance to the RG
    implicit none
    real(fp), intent(in) :: r
    real(fp)             :: r_cond
    r_cond = (R0**2 + r**2)**0.5
  end function physical_distance
  subroutine set_up_photons()
    ! function that defines photon field from magnitude values
    use utils, only : &
         &read_data_files,&
         &spacing=>linspace,&
         &interpol=>simple_approximation
    implicit none
    integer :: l,u
    !                                                                B,       V,       R,       I bands
    real(fp), dimension(4), parameter :: F0                =(/ 4.26e+3, 3.64e+3, 3.08e+3, 2.55e+3 /) ! Jy
    real(fp), dimension(4), parameter :: lambda            =(/ 0.44e-4, 0.55e-4, 0.64e-4, 0.79e-4 /) ! cm
    real(fp), dimension(4), parameter :: dlambda_to_lambda =(/ 0.22,    0.16,    0.23,    0.19    /)

    real(fp), dimension(:,:), allocatable :: B,V,R,I

    call read_data_files(B,"observations/RSOph_Bmag.dat",3,optional_max_length=600)
    call read_data_files(V,"observations/RSOph_Vmag.dat",3,optional_max_length=1200)
    call read_data_files(R,"observations/RSOph_Rmag.dat",3,optional_max_length=600)
    call read_data_files(I,"observations/RSOph_Imag.dat",3,optional_max_length=600)

    ! define time sequence
    call spacing(Time,0._fp,T_max)    

    do l=1,N_t
       mag(1,l)=Time(l)
       mag(2,l) = exp(interpol(Time(l)/86400.,B(1,:),log(B(2,:)))) ! needs to be converted to days
       mag(3,l) = exp(interpol(Time(l)/86400.,V(1,:),log(V(2,:))))
       mag(4,l) = exp(interpol(Time(l)/86400.,R(1,:),log(R(2,:))))
       mag(5,l) = exp(interpol(Time(l)/86400.,I(1,:),log(I(2,:))))
    end do
    Fearth(1,:) = mag(1,:)
    do l=2,5
       Fearth(l,:) = 1.51 * 10 ** (3 - 0.4 * mag(l,:)) * F0(l-1) * (5.03e+15 * lambda(l-1))  ! 1/erg/cm^2/s
    end do
  end subroutine set_up_photons
  
  function get_photons(t,r) result(density)
    ! returns the denstiy of target photons
    use utils, only : erg_to_mcc,h=>hPlanck,c=>light_velocity,spacing=>logspace,AU_to_cm,cm_per_pc
    implicit none
    integer, intent(in)           :: t
    real(fp), intent(in), optional :: r
    integer             :: l
    !                                                                B,       V,       R,       I bands
    real(fp), dimension(4), parameter :: F0                =(/ 4.26e+3, 3.64e+3, 3.08e+3, 2.55e+3 /) ! Jy
    real(fp), dimension(4), parameter :: lambda            =(/ 0.44e-4, 0.55e-4, 0.64e-4, 0.79e-4 /) ! cm
    real(fp), dimension(4), parameter :: dlambda_to_lambda =(/ 0.22,    0.16,    0.23,    0.19    /)
    !
    logical, save :: init_photons = .true.
    real(fp), dimension(size(soft_photon))     :: density
    real(fp)                                   :: distance

    if (init_photons) then
       init_photons = .false.
       soft_photon(N_ph-3:N_ph) = h * c / lambda(4:1:-1) * erg_to_mcc
       call spacing(soft_photon(1:N_ph-3),soft_photon(N_ph-3)/1.e2,soft_photon(N_ph-3))
    end if
    if (present(r)) then
       distance = r
    else
       distance = shock_radius(t)
    endif
    density = 0.
    density(N_ph-3:N_ph) = Fearth(5:2:-1,t) * ( D_to_Earth*cm_per_pc/distance)**2 / c / erg_to_mcc
    density(1:N_ph-3) = Fearth(5,t) * ( D_to_Earth*cm_per_pc/distance)**2 / c / erg_to_mcc 
  end function get_photons
  
  function get_Bfield_time(t) result(B)
    !magnetic field in the downstream
    implicit none
    integer, intent(in) :: t
    real(fp)            :: B
    B = 4 * B_s * (R_s / physical_distance(shock_radius(t)))**2! monopole like B-field compressed by a strong shock by a factor of 4
  end function get_Bfield_time

  function get_Bfield(r) result(B)
    !magnetic field in the downstream
    implicit none
    real(fp), intent(in) :: r
    real(fp)            :: B
    B = 4 * B_s * (R_s / r)**2! monopole like B-field compressed by a strong shock by a factor of 4
  end function get_Bfield

  function get_density_time(t) result(rho)
    ! upstream density
    implicit none
    integer, intent(in) :: t
    real(fp)            :: rho
    real(fp)            :: rrg,v_wind
    rrg      = physical_distance(shock_radius(t))
    v_wind   = beta_speed(rrg)
    rho = dot_M / (4 * pi * v_wind * rrg ** 2 ) 
  end function get_density_time

  function get_density(rrg) result(rho)
    ! upstream density
    implicit none
    real(fp), intent(in) :: rrg
    real(fp)             :: rho
    real(fp)             :: v_wind
    v_wind   = beta_speed(rrg)
    rho = dot_M / (4 * pi * v_wind * rrg ** 2 ) 
  end function get_density
  
end module setup
