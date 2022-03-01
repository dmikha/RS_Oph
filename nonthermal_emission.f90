program nonthermal
  use utils, only : fp=>dp, newunit,read_data_files,linspace,logspace,&
       &TeV_to_mcc,erg_to_mcc,h=>hPlanck,c=>light_velocity,cm_per_pc,AU_to_cm,pi,&
       &interpol=>simple_approximation,&
       &integration=>simple_integration
  use myparameters
  use setup
  use ic, only : ic_cross=>interaction_angle_averaged_ic_cross_section_function
  use pp
  use pair, only : tau_pls=>optical_depth_point_like_photon_source_array
  implicit none
  integer, parameter  :: N_g   = 140, N_x = 120, N_th=100
  integer             :: t,ne,ng,u,nth
  real(fp)            :: r,rrg,B,down_rho,theta_an
  real(fp), dimension(15,N_t) :: physical_condition      ! (1,:) time elapsed since explosion
                                                         ! (2,:) distance from explosion;             
                                                         ! (3,:) distance from RG;                    
                                                         ! (4,:) upstream density;                    
                                                         ! (5,:) upstream magnetic field;
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
  real(fp), dimension(N_e)                  :: Ee,dNe,Ep,dNp,sigma_pp
!  real(fp), parameter                       :: Ee_min = 1.e-4 * TeV_to_mcc, Ee_max=10 * TeV_to_mcc
!  real(fp), parameter                       :: Ep_min = 1.e2, Ep_max=1.e6
  real(fp), dimension(N_e,N_t+1)            :: electrons, protons
  real(fp), dimension(N_g)                  :: vhe,pp_vhe
  real(fp), dimension(N_g,N_t)              :: ic_emission,pp_emission,attenuation
  real(fp), dimension(N_g)                  :: optical_depth
  real(fp), dimension(N_x)                  :: xray
  real(fp), dimension(N_x,N_t)              :: xray_emission
  real(fp), dimension(N_e,N_g)              :: IC

  ! defines array for get_photons function
  call set_up_photons()
  ! energy arrays
  call logspace(Ee,Ee_min,Ee_max) ! for particles
  call logspace(vhe,Ee_min/1.e2,Ee_max*0.999) ! for gamma
  call logspace(pp_vhe,Ep_min/1.e2,Ep_max)        ! for pp
  call logspace(xray,1.e-6_fp,1.e3_fp) ! for xray
  call logspace(Ep,Ep_min,Ep_max) ! for particles

  open(newunit(u),file=out_dir//'physical_conditions_unform.dat',form="unformatted")
  read(u) physical_condition
  close(u)

  open(newunit(u),file=out_dir//'protons.dat',form="unformatted")
  read(u) protons
  close(u)
  open(newunit(u),file=out_dir//'electrons.dat',form="unformatted")
  read(u) electrons
  close(u)

  ! attenuation for electron energies
  attenuation(:,:) = 0.
  do t=N_t/4, N_t, N_t/6 !
     r = physical_condition(2,t-1)
     den = get_photons(t-1,1._fp)
     do nth=1,N_th
        theta_an = nth*3.1415/N_th
        optical_depth=tau_pls(vhe(:),soft_photon(:),den(:),initial_angle=theta_an)/r
        attenuation(:,t) = attenuation(:,t) + sin(theta_an) * exp(-optical_depth)
     enddo
  enddo
  attenuation = attenuation * 3.1415/(N_th*2)
  open(newunit(u),file=out_dir//'gammagamma_examples_ee_check.dat')
  do ng=1,N_g
     write(u,*) vhe(ng), attenuation(ng,N_t/4:N_t:N_t/6)
  enddo
  close(u)

  
  ! leptons emisson
  do t=N_t/4, N_t, N_t/6 !
     dNe = electrons(:,t)
     r = physical_condition(2,t-1)
     B = physical_condition(5,t-1)
     den = get_photons(t-1,r)
     do ng=1,N_g
        do ne=1,N_e
           IC(ne,ng) = c * &
                &integration(soft_photon(:),den(:)* ic_cross(Ee(ne),soft_photon(:),vhe(ng)))
        enddo
        ic_emission(ng,t) = integration(Ee(:),dNe(:)*IC(:,ng))
!        print *,ic_emission(ng,t)
     enddo
  enddo
  open(newunit(u),file=out_dir//'ic_examples.dat')
  do ng=1,N_g
     write(u,*) vhe(ng), ic_emission(ng,N_t/4:N_t:N_t/6)
  enddo
  close(u)
  open(newunit(u),file=out_dir//'gammagamma_ic_examples.dat')
  do ng=1,N_g
     write(u,*) vhe(ng), ic_emission(ng,N_t/4:N_t:N_t/6) * attenuation(ng,N_t/4:N_t:N_t/6)
  enddo
  close(u)


  ! attenuation for protons energies
  attenuation(:,:) = 0.
  do t=N_t/4, N_t, N_t/6 !
     r = physical_condition(2,t-1)
     den = get_photons(t-1,1._fp)
     do nth=1,N_th
        theta_an = nth*3.1415/N_th
        optical_depth=tau_pls(1836.1527*pp_vhe(:),soft_photon(:),den(:),initial_angle=theta_an)/r
        attenuation(:,t) = attenuation(:,t) + sin(theta_an) * exp(-optical_depth)
     enddo
  enddo
  attenuation = attenuation * 3.1415/(N_th*2)
  open(newunit(u),file=out_dir//'gammagamma_examples.dat')
  do ng=1,N_g
     write(u,*) 1836.1527*pp_vhe(ng), attenuation(ng,N_t/4:N_t:N_t/6)
  enddo
  close(u)
  

  pp_emission=0
  ! proton emission Kelner
  do t=N_t/4, N_t, N_t/6 !
     dNp = protons(:,t)
     down_rho = 4 * physical_condition(4,t-1)/1.67e-24 ! jump conditions
     do ng=1,N_g
        pp_emission(ng,t) = integration(Ep(:),dNp(:) * proton_proton_gamma(pp_vhe(ng),Ep(:))) &
             &/proton_proton_cooling_time(down_rho)
     enddo
  enddo
  open(newunit(u),file=out_dir//'pp_examples_Kelner.dat')
  do ng=1,N_g
     write(u,*) pp_vhe(ng), pp_emission(ng,N_t/4:N_t:N_t/6)
  enddo
  close(u)
  open(newunit(u),file=out_dir//'gammagamma_pp_examples_Kelner.dat')
  do ng=1,N_g
     write(u,*) pp_vhe(ng), pp_emission(ng,N_t/4:N_t:N_t/6) * attenuation(ng,N_t/4:N_t:N_t/6)
  enddo
  close(u)
  ! proton emission KATV
  do t=N_t/4, N_t, N_t/6 !
     dNp = protons(:,t)
     down_rho = 4 * physical_condition(4,t-1)/1.67e-24 ! jump conditions
     do ng=1,N_g
        do ne=1,N_e
           sigma_pp(ne) = dif_sigma_proton_proton_KATV(Ep(ne),pp_vhe(ng),model_in='g')
        enddo
        pp_emission(ng,t) = integration(Ep(:), dNp(:) * sigma_pp(:))*down_rho*3.e10
     enddo
  enddo

  open(newunit(u),file=out_dir//'pp_examples.dat')
  do ng=1,N_g
     write(u,*) pp_vhe(ng), pp_emission(ng,N_t/4:N_t:N_t/6)
  enddo
  close(u)
  open(newunit(u),file=out_dir//'gammagamma_pp_examples.dat')
  do ng=1,N_g
     write(u,*) pp_vhe(ng), pp_emission(ng,N_t/4:N_t:N_t/6) * attenuation(ng,N_t/4:N_t:N_t/6)
  enddo
  close(u)
end program nonthermal
