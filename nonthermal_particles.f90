program nonthermal
  use utils, only : fp=>dp, newunit,read_data_files,linspace,logspace,&
       &TeV_to_mcc,erg_to_mcc,h=>hPlanck,c=>light_velocity,cm_per_pc,AU_to_cm,pi,&
       &interpol=>simple_approximation,&
       &integration=>simple_integration,&
       &mls=>make_latex_string
  use myparameters
  use setup
  use ic, only : ic_losses=>ic_losses_single_photon!(electron_energy,photon_energy)
  use synchrotron, only : syn_losses=>synchrotron_losses_function_iso

  use particles, only : part_evol=>distribution_evolution
  implicit none
  integer, parameter  :: N_l = 300, N_int=100
  integer             :: t,l,u,ne,ng
  real(fp)            :: v,vej,r,rrg,rho,B,dEdt_sh,target_energy_density,max_e,max_p,dt,A_norm
  real(fp), dimension(15,N_t) :: physical_condition      ! (1,:) time elapsed since explosion
                                                         ! (2,:) distance from explosion;             
                                                         ! (3,:) distance from RG;                    
                                                         ! (4,:) upstream density;                    
                                                         ! (5,:) (down!)\sout{up}stream magnetic field;
                                                         ! (6,:) shock speed;                         
                                                         ! (7,:) electron maximum energy;                         
                                                         ! (8,:) proton maximum energy;                         
                                                         ! (9,:) target photons energy density
                                                         ! (10,:) L_ph
                                                         ! (11,:) L_sh
                                                         ! (12,:) dE/dt(max)
                                                         ! (13,:) shocked gass mass
                                                         ! (14,:) pressure in the shocked gas
                                                         ! (15,:) ejecta speed
  real(fp), dimension(size(soft_photon))    :: den
  real(fp), dimension(N_e)                  :: Ee,qe,dNe,Ep,qp,dNp
  real(fp), dimension(N_l)                  :: El,losses
  real(fp), dimension(N_ph,N_l)             :: losses_array
  real(fp), dimension(N_e,N_t+1)            :: electrons, protons
  integer                                   :: print_out=0

  ! defines array for get_photons function and TIME
  call set_up_photons()
  ! energy arrays
  call logspace(Ee,Ee_min,Ee_max) ! Lorentz factors for electron
  call logspace(Ep,Ep_min,Ep_max) ! Lorentz factors for protons
  call logspace(El,1._fp,1.e10_fp) ! for losses (electrons only)

  ! initial condtions for particle densities
  dNe = 0.
  dNp = 0.
  ! big array to save particle distributions at different times
  electrons = 0
  protons = 0
  ! the first columns -- particle Lorentz factors
  electrons(:,1) = Ee
  protons(:,1) = Ep
  max_e = 1
  max_p = 1

  ! definining physical conditions
  do t=1,N_t-1
     ! print out of each 500th step
     if ( t / 500 == t / 500. ) print *,t,"out of ",N_t-1," run shock"
     ! function that returns shock radius, speed, and ejecta speed
     ! r   --radius of the shock
     ! v   --speed of the shock
     ! vej --speed of the ejecta
     r = shock_radius(t, vsh=v,vej=vej)

     ! distance to the RG star
     rrg = physical_distance(r)
     
     ! flux of kinetic energy to the shocked wind
     dEdt_sh = 2 * pi * r**2 * get_density( rrg ) * v ** 3

     ! energy losses
     !magnetic field strength
     B = get_Bfield( rrg )

     ! big array to save all relevant physical information, see above for the description
     physical_condition(1,t)  = Time(t+1)/86400 ! days
     physical_condition(2,t)  = r 
     physical_condition(3,t)  = rrg
     physical_condition(4,t)  = get_density(rrg)! 
     physical_condition(5,t)  = B
     physical_condition(6,t)  = v
     physical_condition(11,t) = dEdt_sh
     physical_condition(13,t) =  12.6 * r**2 * v * physical_condition(4,t)
     physical_condition(14,t) = 0.75 * v**2 * physical_condition(4,t)
     physical_condition(15,t) = vej!
     if (t > 1) then
        physical_condition(13,t) = physical_condition(13,t) + physical_condition(13,t-1)
     endif
  end do
  
  
  ! computing non-thermal particles
  do t=1,N_t-1
     ! print out of each 500th step
     if ( t / 500 == t / 500. ) print *,t,"out of ",N_t-1," run for particles"

     r   = physical_condition(2,t)  
     rrg = physical_condition(3,t)
     rho = physical_condition(4,t)
     B   = physical_condition(5,t)
     v   = physical_condition(6,t) 
     dEdt_sh = physical_condition(11,t)

     !ic losses
     den = get_photons(t,r)
     target_energy_density = integration(soft_photon(:),soft_photon(:)*den(:))
     do l=1,N_l
        losses_array(:,l) = ic_losses(El(l),soft_photon(:))
        losses(l) = integration(soft_photon(:),losses_array(:,l)*den(:))
     end do
     ! adding synchrotron losses
     losses = losses + syn_losses(El,B)
     if (t> 1) then
        max_p = max_p+ B * v**2 /(eta_p * 3.e10) * &
             &(Time(t) - Time(t-1)) * 3.195154814716398e-07 ! 3.195154814716398e-07 = e/m_pc^2
     end if
     if ( Bell_2013_protons ) then
        max_p = min(max_p,5.328944623944474e-19 * eta_esc * v**2 * sqrt(dot_M/beta_speed(rrg)))
     endif
     !max energy electrons -- with syn and IC losses
     if (t> 1) then
        dt = (Time(t) - Time(t-1)) / N_int
        do l=1,N_int
           max_e = max_e + B * v**2 / ( eta_e * 3.e10) * dt *0.0005866792055096207 &
                & -exp(interpol(log(max_e),log(El),log(-losses))) * dt
           ! 0.0005866792055096207 = (e/m_ec^2).cgs
           if (max_e < 1) max_e = 1.
        enddo
     endif
     !particle injection. Protons
     if (t > 1) then
        if (max_p > Ep_min) then
           A_norm =  (2-alpha_p) / 0.0015 * kappa_p * dEdt_sh/(max_p**(2-alpha_p) - Ep_min**(2-alpha_p))
           if (print_out == 0 ) then
              print_out = 1
              print *,"A_norm=", A_norm
           endif
           dNp = dNp + Ep**(-alpha_p)*exp(-(Ep/max_p)**beta_p)*(Time(t) - Time(t-1)) * A_norm

        endif
     endif
     protons(:,t+1) = dNp
     !particle injection+cooling. Electron
     ! t>3 is set to avoid epochs with too strong IC losses, just a quick-and-dirty fix
     if (t > 3) then
        if (max_e > Ee_min) then
           dNe = dNe + Ee**(-alpha_e)*exp(-(Ee/max_e)**beta_e)*(Time(t) - Time(t-1)) *&
                &(2-alpha_e)/8.2e-7 * kappa_e *kappa_c(Time(t)) *dEdt_sh/(max_e**(2-alpha_e) - Ee_min**(2-alpha_e))
        endif
        call part_evol(Ee,dNe,El,losses,(Time(t) - Time(t-1)),electrons(:,t+1))
        dNe = electrons(:,t+1)
     endif
     
     physical_condition(7,t) = max_e
     physical_condition(8,t) = max_p
     physical_condition(9,t) = target_energy_density / erg_to_mcc
     physical_condition(10,t) = target_energy_density / erg_to_mcc * 4*pi * r**2 * c ! photon luminocity 
     physical_condition(12,t) = -exp(interpol(log(max_e),log(El),log(-losses)))

  enddo

  open(newunit(u),file=out_dir//'physical_conditions.dat')
  do l=1, N_t - 1
     write(u,*) physical_condition(:,l)
  enddo
  close(u)
  
  open(newunit(u),file=out_dir//'physical_conditions_unform.dat',form="unformatted")
  write(u) physical_condition
  close(u)
  
  open(newunit(u),file=out_dir//'optical_check.dat')
  do l=1, N_t 
     write(u,*) Fearth(:,l)
  enddo
  close(u)

  open(newunit(u),file=out_dir//'protons.dat',form="unformatted")
  write(u) protons
  close(u)
  open(newunit(u),file=out_dir//'electrons.dat',form="unformatted")
  write(u) electrons
  close(u)
  open(newunit(u),file=out_dir//'electrons_examples.dat')
  do l=1,N_e
     write(u,*) electrons(l,1), electrons(l,N_t/4:N_t:N_t/6)
  enddo
  close(u)
  open(newunit(u),file=out_dir//'protons_examples.dat')
  do l=1,N_e
     write(u,*) protons(l,1), protons(l,N_t/4:N_t:N_t/6)
  enddo
  close(u)


  ! write the table for the LaTeX file
  block
    character(len=20) &
         &s_alpha_e,&
         &s_alpha_p,&
         &s_eta_e,&
         &s_eta_p,&
         &s_beta_e,&
         &s_beta_p,&
         &s_eta_esc,&
         &s_kappa_e,&
         &s_kappa_p,&
         &s_emin_e,&
         &s_emin_p,&
         &s_B0,&
         &s_R0,&
         &s_Mdot,&
         &s_Rorb,&
         &s_distance,&
         &s_vej,&
         &s_mej,ss_mej
    write(s_alpha_e,'(F3.1)') alpha_e
    write(s_alpha_p,'(F3.1)') alpha_p
    write(s_eta_e,'(I3)') nint(eta_e/6.28)!3185307179586
    write(s_eta_p,'(I3)') nint(eta_p/6.28)
    write(s_beta_e,'(F3.1)') beta_e
    write(s_beta_p,'(F3.1)') beta_p
    write(s_eta_esc,'(F6.4)') eta_esc
    write(s_kappa_e,'(I2)') nint(kappa_e*100)
    write(s_kappa_p,'(I2)') nint(kappa_p*100)
    write(s_emin_e,'(F3.1)') Ee_min/1.e4
    write(s_emin_p,'(F2.0)') Ep_min
    write(s_B0,'(F2.0)') B_s
    write(s_R0,'(F4.2)') (R_s/AU_to_cm)
    write(s_Mdot,'(F3.1)') (dot_M/v_w/1.e3*1.e2/1.e11)
    write(s_Rorb,'(F4.2)') (R0/AU_to_cm)
    write(s_distance,'(F3.1)') (D_to_Earth/1.e3)
    write(s_vej,'(I1)') nint(v_ej0 / 1.e8)
    write(s_mej,'(F3.1)') (m_ej/2.e26)
    ss_mej = mls(m_ej/2.e33,n=0)
    open(newunit(u), file=out_dir//"table_input.tex", status="replace")
    !inserting values
    write(u,'(A)') "Acceleration slope electrons& \(\alpha_e\)  & -- &\("//trim(s_alpha_e)//"\)  \\"
    write(u,'(A)') "Acceleration slope protons& \(\alpha_p\)  & \("//trim(s_alpha_p)//"\) &--  \\"
    write(u,'(A)') "Cutoff exponent electrons& \(\beta_e\)  & -- &\("//trim(s_beta_e)//"\)  \\"
    write(u,'(A)') "Cutoff exponent protons& \(\beta_p\)  & \("//trim(s_beta_p)//"\) &--  \\"
    write(u,'(A)') "Fraction of energy in electrons&\(\kappa_e\)  & \(0\) &\("//trim(s_kappa_e)//"\%\)  \\"
    write(u,'(A)') "Fraction of energy in protons&\(\kappa_p\)  & \("//trim(s_kappa_p)//"\%\) &\(0\)  \\"
    write(u,'(A)') "Acceleration efficiency of electrons&\(\eta_e\)  & -- &\("//trim(s_eta_e)//"\pi\)  \\"
    write(u,'(A)') "Acceleration efficiency of protons&\(\eta_p\)  & \("//trim(s_eta_p)//"\pi\) &--  \\"
    write(u,'(A)') "Escape efficiency&\(\eta_{\rm esc}\)  & \("//trim(mls(eta_esc,e=.True.))//"\) &--  \\"
    write(u,'(A)') "Electron low energy cutoff&\(E_{\mathrm{min}}\)&--&\("//trim(mls(Ee_min,n=0,e=.True.))//"m_ec^2\)\\"
    write(u,'(A)') "Proton low energy cutoff&\(E_{\mathrm{min}}\)&\("//trim(s_emin_p)//"m_pc^2\)&--\\"
    write(u,'(A)') "\hline"
    write(u,'(A)') "RG surface magnetic field& \(B_*,\rm\, G\) & \("//trim(s_B0)//"\) & \("//trim(s_B0)//"\)\\"
    write(u,'(A)') "RG radius, AU&\(R_{*},\rm\, AU\) & \("//trim(s_R0)//"\) & \("//trim(s_R0)//"\)\\"
    write(u,'(A)') "RG mass-loss rate &\(\dot{M}/v_{\rm w}, \rm\,\kg\cm^{-1}\) & \("//trim(s_Mdot)//&
         &"\times10^{\cgsi{11}{12}}\) & \("//trim(s_Mdot)//"\times10^{\cgsi{11}{12}}\)\\"
    write(u,'(A)') "WD orbit radius&\(d,\rm\, AU\) & \("//trim(s_Rorb)//"\) & \("//trim(s_Rorb)//"\)\\"
    write(u,'(A)') "Distance from Earth&\(\rm kpc\) & \("//trim(s_distance)//"\) & \("//trim(s_distance)//"\)\\"
    write(u,'(A)') "\hline"
    write(u,'(A)') "Ejecta initial speed & \(v_{\mathrm{ej},0},~ {\rm km\,s}^{-1}\) &\("//trim(s_vej)//&
         &"\,000\)&\("//trim(s_vej)//"\,000\)\\"
     write(u,'(A)') "Ejecta mass & \(m_{\mathrm{ej}}\) & \("//trim(ss_mej)//"M_\odot\) &\("//trim(ss_mej)//&
         &"M_\odot\)\\"
    close(u)  
  end block
  
contains
  function kappa_c(time) result(k)
    ! time dependence of the injection, k=1 means steady injection
    implicit none
    real(fp), intent(in) :: time
    real(fp)             :: k
    !k = 30 / (1 + exp(((time/8.64e4)-1)**2/4))
    !k =  1 + 40 / (1 + exp(((time/8.64e4)-1)**2/2))
    k = 1
  end function kappa_c
  
end program nonthermal
