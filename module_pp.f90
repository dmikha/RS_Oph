module pp
  use utils, only : fp=>dp
  implicit none
  real(fp), parameter :: T_th        =  0.2981011622634539 ! kinetic threshold energy (after Eq. 1 in Kafexhiu et al 2014)     
  real(fp), parameter :: mpcc_to_GeV = 0.9382720881604906
  real(fp), parameter :: m_pi        = 0.14385656538566072 ! mpcc units
  
  interface pion_multiplicity
     module procedure  pion_multiplicity_s, pion_multiplicity_v
  end interface pion_multiplicity
  public dif_sigma_proton_proton_KATV
  public E_gamma_by_pp
  public sigma_pion
  public pion_multiplicity
  public proton_proton_gamma!(Eg,Ep)
  public proton_proton_cooling_time!(Ep,n)
contains
  function pp_cross_section(E_p) result(sigma)
    ! E_p in mpcc units
    ! sigma in cm^2
    implicit none
    real(fp), intent(in)    :: E_p
    real(fp)                :: sigma
    real(fp)                :: T_p, L
    T_p = E_p - 1.
    L = log( T_p / T_th )
    sigma = 0.
    if ( T_p > T_th) then
       sigma = 1.e-27 * ( 30.7 - 0.96 * L + 0.18 * L ** 2 ) *&
            &( 1 - (T_th/T_p) ** 1.9 )**3
    endif
  end function pp_cross_section

  function pp_pi_sigma(E_p) result(sigma)
    !E_p in mpcc units
    !sigma in cm^2
    !Eq.2 from KATV 2014
    implicit none
    real(fp), intent(in)      :: E_p
    real(fp)                  :: sigma
    ! some constant defined in Kafexhiu 2014 around Eq.4 
    !real(fp), parameter       :: gamma,K,Gamma_res,M_res,m_pi 
    real(fp), parameter       :: M_res     = 1.26647697932664350 ! mpcc units
    real(fp), parameter       :: Gamma_res = 0.24129461257220572 ! mpcc units
    real(fp), parameter       :: gamma     = (M_res**2 * (M_res**2 + Gamma_res**2))**0.5
    real(fp), parameter       :: K         = 0.9003163161571062 * M_res * Gamma_res * gamma * &
         &(M_res**2 + gamma)**(-0.5)


    real(fp)                  :: s,s05,f_BW,eta
    sigma = 0.
    if ( E_p - 1 > T_th) then
       s     = 2 * (E_p + 1) ! T_p=E_p-1 is kinetic energy
       s05   = s**0.5
       eta   = ((s - m_pi**2 -4)**2 - 16 * m_pi**2)**0.5 / (2 * m_pi * s05)
       f_BW  = K / ( ((s05 - 1)**2 - M_res**2)**2 + M_res**2 * Gamma_res**2)
       sigma = 7.66e-30 * eta**1.95 * (1+eta+eta**5)*f_BW**1.86
    endif
  end function pp_pi_sigma
  function pp_pipi_sigma(E_p) result(sigma)
    ! approximation for proton KINETIC energy between 0.56 and 2 GeV
    !E_p in mpcc units
    !sigma in cm^2
    !Eq.5 from KATV 2014
    implicit none
    real(fp), intent(in)      :: E_p
    real(fp)                  :: sigma
    real(fp)                  :: T_p
    T_p      = (E_p - 1) * 0.9382720881604906 ! here in GeV units (see Fig. 3 caption in KATV) 
    sigma = 0
    ! approximation range
    if ( T_p > 2 ) then
       sigma = 0.
    elseif (T_p > 0.56 ) then
       sigma = 5.7e-27 / ( 1+exp(-9.3*(T_p - 1.4)))
    endif
    return
  end function pp_pipi_sigma


  function pion_multiplicity_s(E_p,model_in,R) result(n)
    !E_p in mpcc units
    !model_in = 'g' for Geant, 's' for Sibyll, 'p' for Pythia, 'q' for QGSJET
    !n is averaged number of pions
    !Eq. 7 from KATV 2014
    implicit none
    real(fp), intent(in)            :: E_p
    integer,  intent(out), optional :: R
    character(1),  intent(in), optional :: model_in
    real(fp)                       :: n
    real(fp), dimension(5)         :: a
    character(1)                   :: model
    real(fp)                       :: T_p, T_pGeV,xi,Q

    T_p    = E_p -1
    T_pGeV = T_p * 0.9382720881604906
    Q      = T_p - T_th ! mpcc
    xi     = T_p - 3.197366774366684 ! (T_p - 3GeV)/ mpcc definition befor Eq.7 KATV
    ! if (present(model_in)) then
    !    model = model_in
    ! else
    !    model = 'g'
    ! endif
    
    model = 'g'
    if ( present(model_in) ) then
       if (model_in == 's' .or. model_in == 'p' .or. model_in == 'q' ) model = model_in
    endif

    n = 0
    
    if ( .not. ( model == 'g' .or. model == 'p' .or. model == 's' .or. model == 'q') ) then
       !stop 'possible models are GEANT4 (g), PYTHIA8 (p), SIBYLL2.1 (s), or QGSJET (q)'
       if (present(R)) R = 0
       return 
    end if

    if ( model == 'g' .and. T_pGeV < 1 ) then
       !stop 'for GEANT model energy should be above 5GeV'
       if (present(R)) R = 0
       return 
    endif
    if ( model == 'p' .and. T_pGeV < 50 ) then
       !stop 'for PYTHIA model energy should be above 50GeV'
       if (present(R)) R = 0
       return 
    endif
    if ( model == 's' .and. T_pGeV < 100 ) then
       !stop 'for SIBYLL model energy should be above 100GeV'
       if (present(R)) R = 0
       return 
    endif
    if ( model == 'q' .and. T_pGeV < 100 ) then
       !stop 'for QGSJET model energy should be above 100GeV'
       if (present(R)) R = 0
       return 
    endif

    a = (/ 0.728, 0.5960, 0.4910, 0.2503, 0.117 /)
    select case (model)
    case ('g')
       a = (/ 0.728, 0.5960, 0.4910, 0.2503, 0.117 /)
    case ('p')
       a = (/ 0.652, 1.6e-3, 0.4880, 0.1928, 0.483 /)
    case ('s')
       a = (/ 5.436, 0.2540, 7.2e-2, 7.5e-2, 0.166 /)
    case ('q')
       a = (/ 0.908, 9.0e-4, 6.0890, 0.1760, 0.448 /)
    end select

    if (T_pGeV < 5) then
       n = -6.e-3 + 0.237 * Q - 0.023 * Q**2 ! Eq. 6 KATV
    else
       n = a(1) * xi**a(4) * (1+exp(-a(2)*xi**a(5)))*(1-exp(-a(3)*xi**0.25))
    endif
    if (present(R)) R = 1
  end function pion_multiplicity_s
  function pion_multiplicity_v(E_p,model_in,R) result(n)
    ! E_p in mpcc
    !model_in = 'g' for Geant, 's' for Sibyll, 'p' for Pythia, 'q' for QGSJET
    !R integer with returns 1 if energy/model are consistent and 0 if not
    !Eq. 7 from KATV 2014
    implicit none
    real(fp), intent(in)            :: E_p(:)
    integer,  intent(out), optional :: R(:)
    character(1),  intent(in), optional :: model_in
    real(fp)                       :: n(size(E_p))
    character(1)                   :: model
    integer    :: i
    model = 'g'
    if ( present(model_in) ) then
       if (model_in == 's' .or. model_in == 'p' .or. model_in == 'q' ) model = model_in
    endif

    do i=1,size(E_p)
       if (present(R)) then 
          n(i) = pion_multiplicity_s(E_p(i),model,R(i))
       else
          n(i) = pion_multiplicity_s(E_p(i),model)
       endif
    end do
    
  end function pion_multiplicity_v
  

  function sigma_pion(E_p, model_in) result(sigma)
    !E_p in mpcc units
    !model_in = 'g' for Geant, 's' for Sibyll, 'p' for Pythia, 'q' for QGSJET
    !sigma in cm^2
    !section 4 from KATV 2014
    implicit none
    real(fp), intent(in)               :: E_p
    character(1), intent(in), optional :: model_in
    real(fp)                           :: sigma
    real(fp)                           :: T_p, T_transition
    character(1)                       :: model
    integer                            :: R
    model = 'g'
    if ( present(model_in) ) then
       if (model_in == 's' .or. model_in == 'p' .or. model_in == 'q' ) model = model_in
    endif
    T_transition = 106578.89247888947
    select case (model)
    case ('g')
       T_transition = 106578.89247888947
       model = 's' ! since the GEANT model works for E<1.e5GeV (see KATV, section 4, first paragrath), switch to SIBYLL (optional, but this model is between Pythia and QGSJET)
    case ('p')
       T_transition = 53.28944623944473
    case ('s')
       T_transition = 106.57889247888946
    case ('q')
       T_transition = 106.57889247888946
    end select
    T_p = E_p - 1.
    sigma = 0
    if ( T_p <= T_th ) return
    if ( T_p <= 2.131577849577789 ) then
       sigma = pp_pi_sigma(E_p)
       sigma = sigma + pp_pipi_sigma(E_p)
       return
    end if
    if ( T_p < T_transition) then
       sigma = pp_cross_section(E_p) * pion_multiplicity_s(E_p,model_in='g')
       return
    end if
    sigma = pp_cross_section(E_p) * pion_multiplicity_s(E_p,model_in=model)
    return
  end function sigma_pion

  function E_pi_max(E_p) result(EpiLABmax)
    ! E_p mpcc units
    ! Emax mpcc units
    ! Eq. 10 from KATV 2014
    implicit none
    real(fp), intent(in)               :: E_p
    real(fp)                           :: EpiLABmax
    real(fp)                           :: s,s05
    real(fp)                           :: EpiCM,PpiCM,gammaCM,betaCM!,gammaLAB,betaLAB

    s       = 2 * (E_p + 1) ! T_p=E_p-1 is kinetic energy
    s05     = s**0.5
    EpiCM   = (s - 4 + m_pi**2)/(2*s05)
    PpiCM   = sqrt(EpiCM**2 - m_pi**2)
    gammaCM = (E_p + 1) / s05
    betaCM  = sqrt(1 - gammaCM**(-2))
    EpiLABmax = gammaCM * ( EpiCM + PpiCM * betaCM)
    return  
  end function E_pi_max

  function E_gamma_by_pp(E_p,mm) result(Egamma)
    ! E_p mpcc units
    ! Emax mpcc units
    ! Eq. 10 from KATV 2014
    implicit none
    real(fp), intent(in)               :: E_p
    character(3), intent(in), optional :: mm
    real(fp)                           :: Egamma
    character(3)                       :: m
    real(fp)                           :: s,s05
!    real(fp)                           :: EpiCM,PpiCM,EpiLABmax,gammaCM,betaCM,gammaLAB,betaLAB
    real(fp)                           :: EpiLABmax,gammaLAB,betaLAB
    m='max'
    if (present(mm)) then
       if (mm == 'min') m = mm
    endif

    !s       = 2 * (E_p + 1) ! T_p=E_p-1 is kinetic energy
    !s05     = s**0.5
    !EpiCM   = (s - 4 + m_pi**2)/(2*s05)
    !PpiCM   = sqrt(EpiCM**2 - m_pi**2)
    !gammaCM = (E_p + 1) / s05
    !betaCM  = sqrt(1 - gammaCM**(-2))
    EpiLABmax = E_pi_max(E_p) !gammaCM * ( EpiCM + PpiCM * betaCM)

    gammaLAB  = EpiLABmax / m_pi
    betaLAB   = sqrt(1 - gammaLAB**(-2))
    if ( m == 'max' ) then
       Egamma = 0.5 * EpiLABmax * ( 1 + betaLAB)
    else
       Egamma = 0.5 * EpiLABmax * ( 1 - betaLAB)
       if (Egamma <= 0 ) Egamma = 0.25 * EpiLABmax * gammaLAB**(-2)
    endif
    return  
  end function E_gamma_by_pp
  
  function kappa_KATV(E_p) result(kappa) ! Eq.14 from Kafexhiu et al 2014
    ! E_p in mpcc units
    implicit none
    real(fp), intent(in) :: E_p
    real(fp)             :: kappa
    kappa = 3.29 - 0.2 * (E_p-1.)**(-1.5) ! T_p is kinetic energy
  end function kappa_KATV
  
  function mu_KATV(E_p) result(mu)
    ! E_p in mpcc units
    ! Eq.14 from KATV 2014    
    implicit none
    real(fp), intent(in) :: E_p
    real(fp)             :: mu
    real(fp)             :: q
    q  = E_p - 2.0657889247888948 ! (E - mpcc -  GeV )/mpcc, note that T_p is kinetic energy
    mu =  1.25 * q ** 1.25 * exp(-1.25*q)
  end function mu_KATV

  function A_max_KATV(E_p,model_in) result(A_max)
    ! E_p in mpcc units
    ! Eq.12 from KATV 2014    
    implicit none
    real(fp), intent(in)               :: E_p
    character(1), intent(in), optional :: model_in
    real(fp)                           :: A_max
    character(1)                       :: model
    real(fp)                           :: T_p, b0, kappa
    real(fp), dimension(3)             :: b
    model = 'g'
    if ( present(model_in) ) then
       if (model_in == 's' .or. model_in == 'p' .or. model_in == 'q' ) model = model_in
    endif
    T_p = E_p - 1
    select case (model)
    case ('s')
       b0 = 0
       kappa = 1.
       b  = (/ 10.77, 0.4120, 1.264e-2 /)
    case ('p')
       b0 = 0
       kappa = 1.
       b  = (/  9.06, 0.3795, 1.105e-2 /)
    case ('q')
       b0 = 0
       kappa = 1.
       b  = (/ 13.16, 0.4419, 1.439e-2 /)
    case default
       if ( T_p < 1.0657889247888945) then ! <1 GeV Eq. 14 KAVT 2014
          b0    = 5.9
          kappa = 3.29 - 0.2 * T_p ** (-1.5) ! Eq. 14 KAVT 2014
       elseif ( T_p < 5.328944623944474 ) then !5GeV
          b0 = 0
          kappa = 1.
          b = (/  9.53, 0.52, 5.4e-2 /)
       else
          b0 = 0
          kappa = 1.
          b = (/  9.13, 0.35, 9.7e-3 /)
       endif
    end select
    A_max = sigma_pion(E_p)
    if ( T_p < 1.0657889247888945) then
       A_max = A_max * b0  / E_pi_max(E_p)
    else
       A_max = A_max * b(1) * T_p **(-b(2)) * exp(b(3)*log(T_p)**2)
    endif
  end function A_max_KATV


  function F_KATV(E_p,Egamma,model_in) result(F)
    ! E_p in mpcc units
    ! Eq.12 from KATV 2014    
    implicit none
    real(fp), intent(in)               :: E_p,Egamma
    character(1), intent(in), optional :: model_in
    real(fp)                           :: F
    character(1)                       :: model
    real(fp)                           :: T_p
    real(fp)                           :: lambda, alpha, beta, gamma, Ymax,Y, X, C, Eg_max
    model = 'g'
    if ( present(model_in) ) then
       if (model_in == 's' .or. model_in == 'p' .or. model_in == 'q' ) model = model_in
    endif
    T_p = E_p - 1
    select case (model)
    case ('s')
       lambda = 3.55
       alpha  = 0.5
       beta   = 3.6
       gamma  = 1
    case ('p')
       lambda = 3.50
       alpha  = 0.5
       beta   = 4.0
       gamma  = 1
    case ('q')
       lambda = 3.55
       alpha  = 0.5
       beta   = 4.5
       gamma  = 1
    case default
       if ( T_p < 1.0657889247888945) then ! <1 GeV Eq. 14 KAVT 2014
          lambda = 1.
          alpha  = 1.
          beta   = 3.29 - 0.2 * T_p**(-1.5)
          gamma  = 0.
       elseif ( T_p <  4.263155699155578 ) then !4GeV KAVT Table V
          lambda = 3.
          alpha  = 1.
          beta   = 2.45 + mu_KATV(E_p)
          gamma  = 1.45 + mu_KATV(E_p)
       elseif ( T_p <  21.315778495777895 ) then !20GeV KAVT Table V
          lambda = 3.
          alpha  = 1.
          beta   = 4.95 + 1.5 * mu_KATV(E_p)
          gamma  = 1.50 + mu_KATV(E_p)
       elseif ( T_p <  106.57889247888946 ) then !100GeV KAVT Table V
          lambda = 3.
          alpha  = 0.5
          beta   = 4.2
          gamma  = 1.
       else
          lambda = 3.
          alpha  = 0.5
          beta   = 4.9
          gamma  = 1.
       endif
    end select
    Eg_max = E_gamma_by_pp(E_p,'max') 
    Ymax = Eg_max  + m_pi**2 / (4*Eg_max)
    C = lambda * m_pi / Ymax
    Y = Egamma + m_pi**2 / (4*Egamma)
    X = (Y - m_pi) / (Ymax - m_pi)
    F = (1 - X**alpha)**beta / ( 1 + X/C)**gamma
    return
  end function F_KATV
  
  function dif_sigma_proton_proton_KATV(E_p,Egamma,model_in) result(sigma)
    ! E_p, Egamma in mpcc units
    ! Eq.8 from KATV 2014    
    implicit none
    real(fp), intent(in)               :: E_p,Egamma
    character(1), intent(in), optional :: model_in
    real(fp)                           :: sigma
    character(1)                       :: model
    model = 'g'
    if ( present(model_in) ) then
       if (model_in == 's' .or. model_in == 'p' .or. model_in == 'q' ) model = model_in
    endif
    sigma = 0.
    if ( E_p - 1 < T_th .or. Egamma >= E_gamma_by_pp(E_p,'max') .or. Egamma <= E_gamma_by_pp(E_p,'min') ) return
    sigma = A_max_KATV(E_p,model_in=model) * F_KATV(E_p,Egamma,model_in=model)
    return
  end function dif_sigma_proton_proton_KATV


  ! functions below are from Kelner et al.
  elemental function proton_proton_gamma(Eg,Ep) result(F)
    ! Eg,Ep in mp c^2 units
    ! returns density of photons for dEg in mpcc units
    implicit none
    real(fp), intent(in) :: Eg,Ep
    real(fp)             :: F
    real(fp)             :: x,L,B,beta,k,xx
    x = Eg / Ep
    L = log(Ep * 5.109989499961642e-07) ! m_pc^2 to TeV units
    B    =  1.300 + 0.140 * L + 0.011 * L ** 2
    beta = (1.790 + 0.110 * L + 0.008 * L ** 2) ** (-1)
    k    = (0.801 + 0.049 * L + 0.014 * L ** 2) ** (-1)
    xx = x ** beta
    F = 0
    if (x < 1) then
       F = 1 / log(x) - 4 * beta * xx / ( 1 - xx)
       F = F - 4 * k * beta * xx * ( 1 - 2 * xx)  /&
            &(1 + k * xx * ( 1 - xx))
       F = F * B * log(x) / x
       F = F * ( (1 - xx) / (1 + k * xx * (1 - xx) ) ) ** 4
       F = F / Ep ! transition from dN/dx to dN/dEg
    endif
  end function proton_proton_gamma
  elemental function  proton_proton_cooling_time(n,Ep) result(t)
    implicit none
    real(fp), intent(in)            :: n
    real(fp), intent(in),optional   :: Ep
    real(fp)                        :: t
    t = 1.e15 / n
  end function proton_proton_cooling_time

end module pp
