module particles
  use utils, only : dp
  implicit none
  private
  public distribution_evolution
   !usage (energy,density_initial,energy_array,energy_losses,dt,density_final)
contains
  function energy_evolution(energy_init,energy_array,energy_losses,delta_t) result(energy_fin)
    use utils, only : spline_coeficients
    implicit none
    real(dp), intent(in) :: energy_init,energy_array(:),energy_losses(:),delta_t
    real(dp) :: energy_fin
    real(dp) :: a,b,c ! spline coefficients
    integer i
    real(dp) :: d,de_dt,v,w,y1,y2,phi,tmp1,tmp2,tmp3 ! tmp variables
    real(dp), parameter :: epsilon=1.e-3_dp
    real(dp) :: t_reminding,energy_current,energy_next_itteration
    

    energy_current=energy_init
    t_reminding=delta_t


    do while (t_reminding > 0._dp)
       call spline_coeficients(energy_current,energy_array,energy_losses,a,b,c,i)

       de_dt = a + b * energy_current + c * energy_current ** 2

       if ( abs(de_dt*t_reminding/energy_current) < epsilon ) then
          energy_next_itteration=energy_current + de_dt*t_reminding
          t_reminding=0._dp
       else

          if ( abs(c * energy_current * t_reminding) < epsilon / 3. ) then
             if (abs(b*t_reminding) <= epsilon/3. ) then
                !  "c,b=0" checked 12/10/2013
                energy_next_itteration=energy_current+a * t_reminding
                if ( energy_next_itteration > energy_array( i ) .and. &
                     &energy_next_itteration < energy_array( i + 1 ) )  then
                   t_reminding=0._dp 
                else if ( energy_next_itteration < energy_array( i ) ) then
                   energy_next_itteration=&
                        &max(energy_array(i)/(1._dp + epsilon),energy_next_itteration)
                   t_reminding=t_reminding+(energy_current-energy_next_itteration)/a
                else if ( energy_next_itteration > energy_array( i + 1 ) ) then 
                   energy_next_itteration=&
                        &min(energy_array( i + 1 ) * ( 1._dp + epsilon),energy_next_itteration)
                   t_reminding=t_reminding + (energy_current-energy_next_itteration) / a
                endif
             else
                !  "c=,b/=0" checked 12/10/2013

                energy_next_itteration=(a / b + energy_current)*exp( b * t_reminding ) - a / b
                if ( energy_next_itteration > energy_array(i) .and. &
                     &energy_next_itteration < energy_array( i + 1 ) )  then
                   t_reminding=0._dp 
                else if ( energy_next_itteration < energy_array( i ) ) then
                   energy_next_itteration=&
                        &max(energy_array( i )/(1._dp + epsilon),energy_next_itteration)
                   t_reminding=t_reminding-log((a+b*energy_next_itteration)/(a+b*energy_current))/b
                else if ( energy_next_itteration > energy_array(i+1) ) then 
                   energy_next_itteration=&
                        &min(energy_array(i+1)*(1._dp+epsilon),energy_next_itteration)
                   t_reminding=t_reminding-log((a+b*energy_next_itteration)/(a+b*energy_current))/b
                endif
             endif
          else
             d=b ** 2 - 4._dp * a * c
             v=b / ( 2._dp * c )
             if ( sqrt(abs(d)) * t_reminding < epsilon ) then  !  this condition is  OK: see page 102
                !  "d=0" checked 12/10/2013

                if ((energy_current+v)*c*t_reminding >= 1._dp ) then ! solution goes to infty
                   
                   if ( c > 0. ) then 
                      energy_next_itteration=energy_array(i+1)*(1._dp+epsilon)
                   else
                      energy_next_itteration=energy_array(i)*(1._dp-epsilon)
                   endif

                else
                   energy_next_itteration=-v-(energy_current+v-d*t_reminding/(4.*c))/((energy_current+v)*c*t_reminding-1._dp) !this equation may miss a term, introduced on may 20th 2014
                endif
                if ( energy_next_itteration > energy_array(i) .and. &
                     &energy_next_itteration < energy_array(i+1) )  then
                   t_reminding=0._dp 
                else if ( energy_next_itteration < energy_array(i) ) then
                    energy_next_itteration=max(energy_array(i)/(1._dp+epsilon),energy_next_itteration)
                   t_reminding=t_reminding-(1._dp/(energy_current+v)-1._dp/(energy_next_itteration+v))/c
                else if ( energy_next_itteration > energy_array(i+1) ) then 
                   energy_next_itteration=min(energy_array(i+1)*(1._dp+epsilon),energy_next_itteration)
                   t_reminding=t_reminding-(1._dp/(energy_current+v)-1._dp/(energy_next_itteration+v))/c
                endif
                                
             else if (d > 0._dp) then
                w=sqrt(d)/(2._dp*c)
                y1=-v+w
                y2=-v-w
                phi=c*t_reminding*w
                tmp1=exp(phi)+exp(-phi)
                tmp2=(exp(phi)-exp(-phi))/tmp1
                tmp1=(energy_current+v)*tmp2/w
                if ( tmp1 >= 1._dp ) then ! solution goes to infty 
                   if ( c > 0. ) then 
                      energy_next_itteration=energy_array(i+1)*(1._dp+epsilon)
                   else
                      energy_next_itteration=energy_array(i)*(1._dp-epsilon)
                   endif

                else
                   energy_next_itteration=-v+(energy_current+v-w*tmp2)/(1._dp-tmp1)
                endif
                if ( energy_next_itteration > energy_array(i) .and. &
                     &energy_next_itteration < energy_array(i+1) )  then
                   t_reminding=0._dp 
                else if ( energy_next_itteration < energy_array(i) ) then
                   energy_next_itteration=&
                        &max(energy_array(i)/(1._dp+epsilon),energy_next_itteration)
                   t_reminding=t_reminding-&
                        &log((energy_next_itteration-y1)*(energy_current-y2)/(energy_next_itteration-y2)/(energy_current-y1))&
                        &/(2.*w*c)
                else if ( energy_next_itteration > energy_array(i+1) ) then 
                   energy_next_itteration=&
                        &min(energy_array(i+1)*(1._dp+epsilon),energy_next_itteration)
                   t_reminding=t_reminding-&
                        &log((energy_next_itteration-y1)*(energy_current-y2)/(energy_next_itteration-y2)/(energy_current-y1))&
                        &/(2.*w*c)
                endif
             else if (d < 0._dp) then
                w=sqrt(-d)/(2._dp*c)
                phi=c*t_reminding*w
                tmp1=cos(phi)
                tmp2=sin(phi)
                tmp3=(energy_current+v)/w
                tmp3=tmp3/sqrt(1._dp+tmp3**2)
                tmp3=acos(tmp3)
                if ( phi >= tmp3 ) then ! solution goes to infty

                   if ( c > 0. ) then 
                      energy_next_itteration=energy_array(i+1)*(1._dp+epsilon)
                   else
                      energy_next_itteration=energy_array(i)*(1._dp-epsilon)
                   endif
                else
                   energy_next_itteration=-v-w*((energy_current+v)*tmp1+w*tmp2)/((energy_current+v)*tmp2-w*tmp1)
                endif
                if ( energy_next_itteration > energy_array(i) .and. &
                     &energy_next_itteration < energy_array(i+1) )  then
                   t_reminding=0._dp 
                else if ( energy_next_itteration < energy_array(i) ) then
                   energy_next_itteration=&
                        &max(energy_array(i)/(1._dp+epsilon),energy_next_itteration)


                   t_reminding=t_reminding-tmp3/(w*c)

                   tmp3=(energy_next_itteration+v)/w
                   tmp3=tmp3/sqrt(1._dp+tmp3**2)
                   tmp3=acos(tmp3)

                   t_reminding=t_reminding+tmp3/(w*c)

                else if ( energy_next_itteration > energy_array(i+1) ) then 
                   energy_next_itteration=&
                        &min(energy_array(i+1)*(1._dp+epsilon),energy_next_itteration)

                   t_reminding=t_reminding-tmp3/(w*c)

                   tmp3=(energy_next_itteration+v)/w
                   tmp3=tmp3/sqrt(1._dp+tmp3**2)
                   tmp3=acos(tmp3)

                   t_reminding=t_reminding+tmp3/(w*c)
                endif
             end if
          endif

       endif
       energy_current=energy_next_itteration
    enddo
    
    energy_fin=energy_next_itteration
    return
  end function energy_evolution

  subroutine distribution_evolution(energy,density_initial,energy_array,energy_losses,dt,density_final)
    implicit none
    real(dp), intent(in) :: energy(:),density_initial(size(energy)),energy_array(:),energy_losses(:),dt
    real(dp), intent(out) :: density_final(size(energy))
    
    real(dp) :: energy_cooled(size(energy)),density_cooled(size(energy))

    if ( size(energy_losses) /= size(energy_array) ) stop "Array mismatch distribution_evolution"
    if ( size(energy) /= size(density_initial) .or. &
         &size(energy) /= size(density_final) ) stop "Array mismatch distribution_evolution"

    call distribution_cooling(energy,density_initial,energy_array,energy_losses,dt,energy_cooled,density_cooled)
    call distribution_transformation(energy_cooled,density_cooled,energy,density_final)

    return
  end subroutine distribution_evolution
  


  subroutine distribution_cooling(energy,density_initial,energy_array,energy_losses,dt,energy_cooled,density_cooled)
    implicit none
    real(dp), intent(in) :: energy(:),density_initial(:),energy_array(:),energy_losses(:),dt
    real(dp), intent(out) :: energy_cooled(:),density_cooled(:)

! temporaly variables  
    integer :: i
    real(dp), parameter  :: epsilon=1.e-1_dp ! relative size of the interval for which density is computed
    real(dp) :: energy_shifted

    if ( size(energy_losses) /= size(energy_array) ) stop "Array mismatch distribution_cooling"
    if (  size(energy) /= size(density_initial) .or. &
         &size(energy) /= size(density_cooled) .or. &
         &size(energy) /= size(energy_cooled) ) stop "Array mismatch distribution_cooling"

    do i=1,size(energy)
       energy_cooled(i)=energy_evolution(energy(i),energy_array,energy_losses,dt)
    enddo
    do i=1,size(energy)
       if (energy(i)*(1._dp+epsilon) < energy_array(size(energy_array)) ) then
          energy_shifted=energy_evolution(energy(i)*(1._dp+epsilon),energy_array,energy_losses,dt) 
          density_cooled(i)=density_initial(i)*epsilon*energy(i)/(energy_shifted-energy_cooled(i))
       else
          energy_shifted=energy_evolution(energy(i)*(1._dp-epsilon),energy_array,energy_losses,dt)
          density_cooled(i)=density_initial(i)*epsilon*energy(i)/(energy_cooled(i)-energy_shifted)
       endif
    enddo
    return
  end subroutine distribution_cooling


  subroutine distribution_transformation(energy_initial,density_initial,energy,density) ! rescales particle distribution to another energy range
    implicit none
    real(dp), intent(in) :: energy_initial(:),density_initial(:),energy(:)
    real(dp), intent(out) :: density(:)

    integer :: i,j
    real(dp) :: tmp_x1,tmp_x2,tmp_x,tmp_y1,tmp_y2,tmp_y


    if ( size(energy_initial) /= size(density_initial) ) stop "Array mismatch distribution_transformation"
    if ( size(energy) /= size(density) ) stop "Array mismatch distribution_transformation"


    j=1
    do i=1,size(energy)
       do while ( energy(i) >= energy_initial(j+1) .and. j < size(energy_initial) )
          j=j+1
       enddo
       if ( energy(i) < energy_initial(j) .or. energy(i) > energy_initial(size(energy_initial)) ) then
          density(i)=0._dp
       else
          if (density_initial(j) > 0._dp .and. density_initial(j+1) > 0._dp ) then
             tmp_x1=log(energy_initial(j))
             tmp_x2=log(energy_initial(j+1))
             tmp_x=log(energy(i))
             tmp_y1=log(density_initial(j))
             tmp_y2=log(density_initial(j+1))
             
             tmp_y=tmp_y1+(tmp_x-tmp_x1)/(tmp_x2-tmp_x1)*(tmp_y2-tmp_y1) ! linear fit for log-values
             
             density(i)=exp(tmp_y)
          else
             density(i)=0._dp
          endif
       endif
       
    enddo

  end subroutine distribution_transformation

end module particles
