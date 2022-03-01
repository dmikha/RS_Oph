module myparameters
  use utils, only : fp=>dp, AU=>AU_to_cm,pi,TeV_to_mcc
  implicit none
  integer,  public, parameter :: N_t        = 7200 ! number of time steps              20000!
  integer,  public, parameter :: N_e        = 100!, N_l = 300, N_int=100!, N_ll = 30
                                            
  integer,  public, parameter :: N_ph       = 100  ! number of target photon energies
  real(fp), public, parameter :: dot_M      = 2*6.3e18 ! g/s 6.3e18 g/s -> 1.e-7 M_sun / yr
  real(fp), public, parameter :: B_s        = 1. !G From Jonathan B = 0.1 G (1AU/r)^2
  real(fp), public, parameter :: R_s        = 0.35 * AU !cm ! O'Brien et al. (1992)
  real(fp), public, parameter :: v_w        = 2.e6 !cm/s
  real(fp), public, parameter :: R0         = 1.48 * AU ! cm ! radius of the orbit
  real(fp), public, parameter :: D_to_Earth = 1.4e3_fp  ! pc ! distance to the Earth
  real(fp), public, parameter        :: T_max = 5.184e5 ! s 4.32e5=> 5 days ; 4.752e5 => 5.5 days
  real(fp), public, dimension(N_t)   :: Time
  real(fp), public, dimension(5,N_t) :: mag, Fearth ! magnitude and above atmosphere flux
  real(fp), public, dimension(N_ph)  :: soft_photon
  real(fp), public, parameter        :: m_ej   = 2.e26 !26    !g 2.e26g = 1.e-7M_sun
  real(fp), public, parameter        :: v_ej0   = 3.00e8_fp  !cm/s 
  real(fp), public, parameter        :: fill    = 1  !filling factor
  integer,  public                   :: error_occured = 0
  real(fp), public, parameter        :: Ee_min = 1.e2, Ee_max=10 * TeV_to_mcc !1.e-4 * TeV_to_mcc !1.e-4 * TeV_to_mcc
  real(fp), public, parameter        :: Ep_min = 2.e0, Ep_max=1.e6
  character(len=*), public, parameter :: out_dir="dat/"
  ! nonthermal particles
  real(fp), parameter :: kappa_e = 3.e-2_fp !fraction of energy to non-thermal
  real(fp), parameter :: kappa_p = 5.e-1_fp !fraction of energy to non-thermal
  real(fp), parameter :: alpha_e = 2.2_fp !injection slope
  real(fp), parameter :: alpha_p = 2.2_fp !injection slope
  real(fp), parameter :: beta_e  = 0.5_fp !exponent power
  real(fp), parameter :: beta_p  = 0.5_fp  !exponent power
  real(fp), parameter :: eta_e   = 10 * 6.28!5*10 !acceleration efficiency 6.28 -> 2pi (c/v)**2
  real(fp), parameter :: eta_p   = 30 * 6.28!30*30 !acceleration efficiency 6.28 -> 2pi (c/v)**2
  !  real(fp), parameter :: fill    = 1.   !filling factor not used
  logical,  parameter :: Bell_2013_protons = .True. ! Bell 2013 escape limit for protons
  real(fp), parameter :: eta_esc = 1.e-2_fp ! Bell 2013 CR flux constnat
  
end module myparameters
