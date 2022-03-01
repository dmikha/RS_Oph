module utils
  implicit none
  private
  integer, parameter, public:: dp=kind(0.d0)
  integer, parameter, public:: sp=kind(0.)
  !constants
  real(dp), parameter, public:: pi=3.141592653589793
  !physical constants in cgs units
  real(dp), parameter, public:: light_velocity=3.e10, stefan_boltzmann=5.6704e-5, elementary_charge=4.8e-10, &
       &electron_mass=9.10938291e-28, electron_radius=2.8179e-13, sigma_t=6.65245873498e-25, &
       &AU_to_cm=1.49598e+13, K_to_mcc=0.1686e-9,&
       &erg_to_mcc=1.221433e+06,erg_to_eV=6.24150934326018e+11,&
       &TeV_to_mcc=1.956947e+06,GeV_to_mcc=1.956947e+03,&
       &keV_to_mcc=1.956947e-03,eV_to_mcc=1.956947e-06,Hz_to_mcc=8.09329978572e-21,&
       &cm_per_pc=3.085680e+18,cm_per_Mpc=3.085680e+24,&
       &s_per_yr=3.155760e+07,s_per_day=8.64e4,&
       &hbar=1.054e-27,hPlanck=6.62607015e-27,&
       &gravitational_contant=6.673e-08,sun_mass=1.989100e+33


  public &
       &make_latex_string,&
       &logspace,&
       &linspace,&
       &linspace_middle_point,&
       &simple_integration,&
       &pl_integration,&
       &simple_approximation,&
       &spline_coeficients,&
       &newunit,&
       &read_data_files
contains
    function value_index_nonstrict(A,aa) result(na) ! if aa<=A(1) then reterns 1, otherwise na>1
      implicit none
      real(dp), intent(in) :: A(:),aa
      integer :: na
      na=1
      do while ( aa > A(na) .and. na<size(A) )
         na=na+1
      end do
      return
    end function value_index_nonstrict

  function make_latex_string(x,n,e) result(out_string)
    implicit none
    real(dp), intent(in) :: x
    integer, optional    :: n
    logical, optional    :: e
    character(len=30)    :: out_string
    integer              :: nd
    logical              :: ei
    character(len=1)     :: s_nd
    character(len=6)     :: s_print
    ! return zero
    if ( x == 0.) then
       out_string = "0"
       return
    end if
    ! check exponential output
    if (present(e)) then
       ei = e
    else
       ei = .False.
    endif
    if ( .not. ei) then
       !number of digits for the
       if (present(n)) then
          nd = n
       else
          nd = 2
       endif
       write(s_nd,'(I1)') nd
       s_print = "(F5."//s_nd//")"
       if ( x > 0.01 .and. x < 10) then
          write(out_string,s_print) x
          return
       end if
       if ( x > -10 .and. x < -0.01) then
          write(out_string,s_print) x
          return
       end if
       if ( x > -100 .and. x <= -10) then
          write(out_string,'(F5.1)') x
          return
       end if
       if ( x >= 10 .and. x < 100) then
          write(out_string,'(F5.1)') x
          return
       end if
       if ( x > -1000 .and. x <= -100) then
          write(out_string,'(F5.0)') x
          return
       end if
       if ( x >= 100 .and. x < 1000) then
          write(out_string,'(F5.0)') x
          return
       end if
    endif
    block
      character(3)  :: s_coef, s_power, s_sign
      real(dp)      :: coef, logg
      integer       :: power
      ! number of digits
      if (present(n)) then
         nd = n
      else
         nd = 1
      endif
      write(s_nd,'(I1)') nd
      s_print = "(F3."//s_nd//")"
       
      if (x > 0) then
         s_sign = ""
      else
         s_sign = "-"
      endif
      logg = log10(abs(x))
      power = nint(logg)
      coef = abs(x) / (10.** power)
      if ( coef > 0.95 .and. coef < 1.05 ) then
         write(s_power,'(I2)') power
         out_string = trim(s_sign)//"10^{"//trim(s_power)//"}"
      else
         if (power <= 0) power = power - 1
         coef = abs(x) / (10.** power)
         write(s_coef,s_print) coef
         write(s_power,'(I2)') power
         if ( power == 0 ) then
            out_string = trim(s_sign)//trim(s_coef)
         else
            out_string = trim(s_sign)//trim(s_coef)//"\times10^{"//trim(s_power)//"}"
         endif
      endif
    end block
    return
  end function make_latex_string
  subroutine linspace(E,Emin,Emax)
    real(dp) E(:)
    real(dp) Emin,Emax
    integer i
    E=[(Emin+(Emax-Emin)*(i-1._dp)/(size(E)-1._dp), i=1,size(E))]
    return
  end subroutine linspace

  subroutine linspace_middle_point(E,Emin,Emax)
    real(dp) E(:)
    real(dp) Emin,Emax
    integer i
    E=[(Emin+(Emax-Emin)*(i-0.5_dp)/size(E), i=1,size(E))]
    return
  end subroutine linspace_middle_point

  function simple_integration(x,y) result(sum)
    real(dp), intent(in)::x(:),y(:)
    real(dp) sum
    integer k
    if (size(x) /= size(y)) stop "inconsistent size of arrays in simple_integration"
    sum=0._dp
    do k=1,size(x)-1
       sum=sum+(y(k)+y(k+1))/2._dp*(x(k+1)-x(k))
    enddo
    return
  end function simple_integration

    function pl_integration(x,y) result(sum)
      ! function reterns an integral assuming a power-law interpolation between the nods
      real(dp), intent(in)::x(:),y(:)
      real(dp) alpha,coef,sum
      integer k
      if (size(x) /= size(y)) stop "inconsistent size of arrays in pl_integration"
      if (size(x) < 2) stop "inconsistent size of arrays in pl_integration"
      if (minval(x) <= 0) stop "inconsistent values in pl_integration, x"
      if (minval(y) < 0) stop "inconsistent values in pl_integration, y"
      sum=0._dp
      do k=1,size(x)-1
         if (y(k+1) > 0 .and. y(k) > 0) then
            alpha = log(y(k+1)/y(k)) / log(x(k+1)/x(k))
            if (alpha > -1.01 .and. alpha < -0.99 ) alpha = -1
            coef = y(k) / x(k)**alpha
            if (alpha == -1) then
               sum=sum + coef * log(x(k+1)/x(k))
            else
               sum=sum + coef / (alpha+1) * (x(k+1)**(alpha+1) - x(k)**(alpha+1))
            endif
         endif
      enddo
      return
    end function pl_integration
  
    function simple_approximation(x,x_array,y_array) result(y)
      implicit none
      real(dp), intent(in) :: x,x_array(:),y_array(:)
      real(dp) :: y
      integer :: i
      if( size(x_array) /= size(y_array) .or. size(x_array) < 2 ) stop "Incorrect size of x,y arrays in simple arrpoximation"
      i=value_index_nonstrict(x_array,x)
      if ( i == 1 ) i=2
      y=y_array(i-1)+(x-x_array(i-1))/(x_array(i)-x_array(i-1))*(y_array(i)-y_array(i-1))
      return
    end function simple_approximation

    integer function newunit(unit) result(n)
    ! returns lowest i/o unit number not in use
    integer, intent(out), optional :: unit
    logical inuse
    integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
    integer, parameter :: nmax=999  ! may be system-dependent
    do n = nmin, nmax
       inquire(unit=n, opened=inuse)
       if (.not. inuse) then
          if (present(unit)) unit=n
          return
       end if
    end do
    stop "newunit ERROR: available unit not found."
  end function newunit

  subroutine read_data_files(array_to_fill,data_file,number_of_raws,&
       &optional_max_length,optional_comment_char)
    ! array_to_fill == double precision array, allocatable
    ! string with the file name
    ! number of raws in the data file
    ! max length of the file, if not given set to 100
    ! comment character, if not given set to #
!    use utils, only : newunit
    implicit none
    real(dp), dimension(:,:), intent(out), allocatable :: array_to_fill
    character(len=*), intent(in) :: data_file
    integer, intent(in) :: number_of_raws
    integer, intent(in), optional :: optional_max_length
    character(1), intent(in), optional :: optional_comment_char
    character(1) :: comment_char
    real(dp), dimension(:,:), allocatable :: read_data
    integer :: max_length
    logical :: was_read=.false.
    integer :: u,stat,number_data_points,row
    character (len=1000) :: line_from_file
    if ( present(optional_max_length)) then
       max_length=optional_max_length
    else
       max_length=100
    endif
    if ( present(optional_comment_char)) then
       comment_char=optional_comment_char
    else
       comment_char="#"
    endif
    allocate(read_data(number_of_raws,max_length))
    print *,"reading data from ",data_file
    open(newunit(u), file=data_file)
    row=0
    do while (row < max_length)
       read (u,'(A)',iostat=stat) line_from_file
       if ( stat > 0 ) then
          stop 'An error occured while reading file'
       elseif ( stat < 0 ) then
          number_data_points = row
          print *, 'EOF of ',data_file,' reached. Found a total of ', number_data_points, 'rows.'
          exit
       endif
       if ( index(line_from_file, comment_char) == 0 ) then
          row=row + 1
          read (line_from_file, *) read_data(:,row)
       endif
    end do

    if ( .not. stat < 0 ) then
       number_data_points = row
       print *, 'EOF of ',data_file,' NOT reached. Read a total of ', number_data_points, 'rows.'
    endif
    close(u)
    allocate(array_to_fill(number_of_raws,number_data_points))
    array_to_fill=read_data(:,:number_data_points)
  end subroutine read_data_files
  subroutine logspace(E,Emin,Emax)
    real(dp) E(:)
    real(dp) Emin,Emax
    integer i
    E=[(Emin*exp(log(Emax/Emin)*(i-1._dp)/(size(E)-1._dp)),i=1,size(E))]
    return
  end subroutine logspace
  function logspace_f(Ar,N,Amin,Amax) result(A_out)
    real(dp), optional :: Ar(:)
    integer, optional :: N
    real(dp) Amin,Amax
    real(dp), dimension(:), allocatable :: A_out
    integer N_size,i
    if ( .not. (present(Ar) .or. present(N)) ) stop "logspace_f error"
    if (present(Ar)) then
       N_size=size(Ar)
    else
       N_size=N
    endif
    allocate(A_out(N_size))
    A_out=[(Amin*exp(log(Amax/Amin)*(i-1._dp)/(N_size-1._dp)),i=1,N_size)]
    return
  end function logspace_f
    subroutine spline_coeficients(x,x_array,y_array,a,b,c,middle_array_index) ! this is not a spline, the overall line can contain breaks, which might be reflected in the subroutine 
      ! calculating distribution evolution

      ! x_array -- array of variable values
      ! y_array -- array of function values
      ! x -- variable value
      ! middle_array_index -- the previous to x node of array
      ! a,b,c quadratic spline coefficients

      implicit none
      real(dp), intent(in) :: x,x_array(:),y_array(:)
      real(dp), intent(out) :: a,b,c !constant, linear, and quadratic spline coefficints
      integer, intent(out) :: middle_array_index !location in the array

      integer n1,n2,n3 ! for quadratic spline

      if( size(x_array) /= size(y_array) .or. size(x_array) < 2 ) stop "Incorrect size of x,y arrays in spline_coeficients"
      if ( x < x_array(1) .or. x>x_array(size(x_array)) ) then
         write(*,*) x,x_array(1)/x,x_array(size(x_array))/x,"from 'spline_coeffients'"
         stop "x is out of range in spline_coeficients"
      endif
      select case (size(x_array))
      case (2)
         if ( x_array(1) == x_array(2) )  stop "Incorrect x-array in spline_coeficients"
         a=(y_array(1)*x_array(2)-y_array(2)*x_array(1))/(x_array(2)-x_array(1))
         b=(y_array(2)-y_array(1))/(x_array(2)-x_array(1))
         c=0._dp
         middle_array_index=1
      case default
         n1=1
         do while (x >= x_array(n1+1) ) ! x should be close to the central point
            n1=n1+1
         enddo

         if (n1==1 .or. &
              &( x > sqrt(x_array(n1)*x_array(n1+1)) .and. n1 < size(x_array)-1 ) ) then 
            n2=n1+1
            n3=n1+2
         else
            n2=n1+1
            n3=n1-1
         endif

         
         a=&
              &y_array(n1)*x_array(n2)*x_array(n3)/(x_array(n1)-x_array(n2))/(x_array(n1)-x_array(n3))+&
              &y_array(n2)*x_array(n1)*x_array(n3)/(x_array(n2)-x_array(n1))/(x_array(n2)-x_array(n3))+&
              &y_array(n3)*x_array(n2)*x_array(n1)/(x_array(n3)-x_array(n2))/(x_array(n3)-x_array(n1))
         b=-1._dp*(&
              &y_array(n1)*(x_array(n2)+x_array(n3))/(x_array(n1)-x_array(n2))/(x_array(n1)-x_array(n3))+&
              &y_array(n2)*(x_array(n1)+x_array(n3))/(x_array(n2)-x_array(n1))/(x_array(n2)-x_array(n3))+&
              &y_array(n3)*(x_array(n2)+x_array(n1))/(x_array(n3)-x_array(n2))/(x_array(n3)-x_array(n1)))

         c=&
              &y_array(n1)/(x_array(n1)-x_array(n2))/(x_array(n1)-x_array(n3))+&
              &y_array(n2)/(x_array(n2)-x_array(n1))/(x_array(n2)-x_array(n3))+&
              &y_array(n3)/(x_array(n3)-x_array(n2))/(x_array(n3)-x_array(n1))

         middle_array_index=n1
      end select
      return
    end subroutine spline_coeficients

 !! Array manipulation

! ! needed for PGI fortran only  requires additionally  -Kieee
!     ! logical function isnan(a)  result(check)
!     !   real(dp), intent(in) ::  a
!     !   logical :: check
!     !   if (a.ne.a) then
!     !      check = .true.
!     !   else
!     !      check = .false.
!     !   endif
!     !   return
!     ! end function isnan
    
 
!     subroutine logspace_middle_point(E,Emin,Emax)
!       real(dp) E(:)
!       real(dp) Emin,Emax
!       integer i
!       E=[(Emin*exp(log(Emax/Emin)*(i-0.5_dp)/size(E)),i=1,size(E))]
!       return
!     end subroutine logspace_middle_point


!     function logspace_function(Emin,Emax,N) result(E)
!       integer, intent(in) :: N
!       real(dp), intent(in) :: Emin,Emax
!       real(dp) ::  E(N)
!       integer i
!       E=[(Emin*exp(log(Emax/Emin)*(i-1._dp)/(N-1._dp)),i=1,N)]
!       return
!     end function logspace_function

!     function logspace_middle_point_function(Emin,Emax,N) result(E)
!       integer, intent(in) :: N
!       real(dp), intent(in) :: Emin,Emax
!       real(dp) ::  E(N)
!       integer i
!       E=[(Emin*exp(log(Emax/Emin)*(i-0.5_dp)/N),i=1,N)]
!       return
!     end function logspace_middle_point_function



!     function  linspace_function(Emin,Emax,N) result(E)
!       integer, intent(in) :: N
!       real(dp), intent(in) ::  Emin,Emax
!       real(dp) :: E(N)
!       integer i
!       E=[(Emin+(Emax-Emin)*(i-1._dp)/(N-1._dp), i=1,N)]
!       return
!     end function linspace_function

!     function Heaviside_function_ar(A) result(B)
!       implicit none
!       real(dp), intent(in) :: A(:)
!       real(dp) :: B(size(A))
!       integer :: i
!       B=0.e0_dp
!       do i=1,size(A)
!          if ( A(i) > 0.e0_dp ) B(i)=1.e0_dp
!       enddo
!       return
!     end function Heaviside_function_ar
!     function Heaviside_function_non_ar(A) result(B)
!       implicit none
!       real(dp), intent(in) :: A
!       real(dp) :: B
!       B=0.e0_dp
!       if ( A > 0.e0_dp ) B=1.e0_dp
!       return
!     end function Heaviside_function_non_ar
    
      
    
!     subroutine array_merge(A,B,C,Range)
!       !subroutine declaration variables
!       real(dp), intent(in) :: A(:,:),B(:,:)
!       real(dp), intent(in out) :: C(:,:)
!       integer, optional :: Range !1 means C(1,:) is defined, 2 means C(1,1) and C(1,-1) are defined the rest is to be log-spaced, others log-spacing other entire interval
!       !internal variables
!       integer Na,Nb,Nc
!       integer :: ir,i,j,k
!       real(dp) :: Emin,Emax,E

!       if (size(A(:,1)) /= 2 .or. size(B(:,1)) /= 2 .or. size(C(:,1)) /= 2) stop "Incorrect size of 'A' or 'B' or 'C' &
!            &in array merge subroutine"

!       Na=size(A(1,:))
!       Nb=size(B(1,:))
!       Nc=size(C(1,:))

!       ! energy interval selection
!       if (present(Range)) then
!          ir=Range
!       else
!          ir=0
!       endif
!       select case (ir)
!       case (1)
!       case (2) !boundary values are defined
!          Emin=C(1,1)  !  
!          Emax=C(1,Nc) ! 
!          call logspace(C(1,:),Emin,Emax)
!       case default 
!          Emin=min(A(1,1),B(1,1))!the lower end
!          Emax=max(A(1,Na),B(1,Nb))!the higher end 
!          call logspace(C(1,:),Emin,Emax)
!       end select
!       ! end of energy interval selection

!       do i=1,Nc
!          E=C(1,i)
!          C(2,i)=0.d0
!          if (E>=A(1,1).and.E<A(1,Na)) then
!             j=1
!             do while (E>A(1,j+1))
!                j=j+1
!             enddo
!             C(2,i)=A(2,j)+(E-A(1,j))/(A(1,j+1)-A(1,j))*(A(2,j+1)-A(2,j))
!          endif

!          if (E>=B(1,1).and.E<B(1,Nb)) then
!             j=1
!             do while (E>B(1,j+1))
!                j=j+1
!             enddo
!             C(2,i)=C(2,i)+B(2,j)+(E-B(1,j))/(B(1,j+1)-B(1,j))*(B(2,j+1)-B(2,j))
!          endif
!       enddo

!       return
!     end subroutine array_merge


    
    
    
!     subroutine chi_squar(y_data,dy_data,y,chi2,norm,x_data,x,log_interpolation)
!       implicit none
!       real(dp), intent(in) :: y_data(:),dy_data(:),y(:)
!       real(dp), intent(inout) :: norm
!       real(dp), intent(out) :: chi2
!       real(dp), intent(in), optional :: x_data(:),x(:)
!       logical, intent(in), optional :: log_interpolation
!       real(dp), dimension(size(y_data)) :: y0,err,y0_data,dy0_data
!       logical :: log_int
!       integer :: k
!       real(dp) :: sum1,sum2

!       if (present(log_interpolation)) then
!          log_int=log_interpolation
!       else
!          log_int=.true.
!       endif

!       if ( present(x) .and. present(x_data) ) then
!          if (size(y_data) /= size(x_data) ) stop "inconsistent size of arrays in chi_squar: y_data and x_data"         
!          if (size(x) /= size(y) ) stop "inconsistent size of arrays in chi_squar: y and x"         
!          if ( log_int ) then
!             !y0 = exp( (/ (polynomial_interpolotation(x_data(k),x,log(y),inter_power=2,error=err(k)),k=1,size(x_data) ) /) )
!             y0 = exp( (/ (simple_approximation(x_data(k),x,log(y)),k=1,size(x_data) ) /) )
!          else
!             !y0 = (/ (polynomial_interpolotation(x_data(k),x,y,inter_power=2,error=err(k)),k=1,size(x_data) ) /) 
!             y0 = (/ (simple_approximation(x_data(k),x,y),k=1,size(x_data) ) /) 
!          endif
! !         print *,y0 / x_data**2,err
!       else
!          if (size(y_data) /= size(y) ) stop "inconsistent size of arrays in chi_squar: y_data and y"         
!          y0=y
!       endif
!       where ( dy_data == 0 )
!          y0 = 0
!          y0_data = 0
!          dy0_data = 1.
!       elsewhere
!          dy0_data = dy_data
!          y0_data = y_data
!       endwhere
!       ! do k=1,size(dy_data)
!       ! if (dy_data(k) /= dy0_data(k) ) then
!       !    print *,"some points excluded"
!       ! end if
!       ! enddo
!       sum1=sum(y0**2 / dy0_data ** 2 )
!       sum2=sum(y0*y0_data / dy0_data ** 2 )
!       if ( norm == 0 ) then
!          norm = sum2 / sum1
!       endif
!       chi2 = sum( (y0_data - y0*norm) ** 2 / dy0_data ** 2 )
!     end subroutine chi_squar
    
!     subroutine powerlaw_fit(x,y,norm,index,x0,x1)
!       ! y=norm * x ** (-index)
!       ! Numerical recipes p. 780
!       implicit none
!       real(dp), intent(in) :: x(:),y(:)
!       real(dp), intent(out) :: norm,index
!       real(dp), intent(in), optional :: x0,x1
!       real(dp) :: e1,e2
!       integer :: i1,i2,n
!       real(dp), dimension(:), allocatable :: e0,d0
!       real(dp) :: S,Sx,Sy,Sxx,Syy,Sxy,Delta
!       if (size(x) /= size(y) .or. size(x) == 1 ) stop "inconsistent size of arrays in powerlaw_fit"
!       if (present(x0)) then
!          i1=value_index(x,x0)
!       else
!          i1=1
!       endif
!       if (present(x1)) then
!          i2=value_index(x,x1)
!       else
!          i2=size(x)
!       end if

!       if ( i1 > i2 ) then
!          i1=i1+i2
!          i2=i1-i2
!          i1=i1-i2
!       end if

!       if (i1 > 1) i1=i1-1
      
!       if (i2 == i1 .and. i2 < size(x) ) i2=i2+1

!       if ( (i2-i1) == 1 ) then
!          index=-log(y(i2)/y(i1))/log(x(i2)/x(i1))
!          norm=y(i2) * x(i2) ** index
!       else
!          n=i2-i1+1
!          allocate(e0(n))
!          allocate(d0(n))
!          e0=log(x(i1:i2))
!          d0=log(y(i1:i2))
!          !       d0=log(energy(i1:i2) ** -1.3 * exp(-energy(i1:i2)/energy(i2)))
!          S=n
!          Sx=sum(e0)
!          Sy=sum(d0)
!          Sxx=sum(e0 ** 2)
!          Syy=sum(d0 ** 2)
!          Sxy=sum(e0 * d0)
!          Delta=S * Sxx - Sx ** 2
!          norm=exp((Sxx * Sy - Sx * Sxy) / Delta)
!          index=(-S * Sxy + Sx * Sy) / Delta
!          !fit(3)=sqrt(See/Delta) !  here is not clear
!          !fit(4)=sqrt(S/Delta)
!       endif
!     end subroutine powerlaw_fit
    
    
!     function array_has_nan_element1D(A) result(check)
!       real(dp), intent(in) :: A(:)
!       logical check
!       integer :: k
!       check=.false.
!       do k=1,size(A)
!          if(isnan(A(k))) check=.true.
!       enddo
!       return
!     end function array_has_nan_element1D

!     function array_has_nan_element2D(A) result(check)
!       real(dp), intent(in) :: A(:,:)
!       logical check
!       integer :: k,i
!       check=.false.
!       do k=1,size(A(1,:))
!          do i=1,size(A(:,1))
!             if(isnan(A(i,k))) check=.true.
!          enddo
!       enddo
!       return
!     end function array_has_nan_element2D

!     function array_has_nan_element3D(A) result(check)
!       real(dp), intent(in) :: A(:,:,:)
!       logical check
!       integer :: k,i,j
!       check=.false.
!       do k=1,size(A(1,1,:))
!          do i=1,size(A(1,:,1))
!             do j=1,size(A(:,1,1))
!                if(isnan(A(j,i,k))) check=.true.
!             enddo
!          enddo
!       enddo
!       return
!     end function array_has_nan_element3D

!     logical function isinf(a) 
!       real(dp), intent(in) ::  a 
!       if ((a*0).ne.0) then 
!          isinf = .true. 
!       else 
!          isinf = .false. 
!       end if
!       return 
!     end function isinf


!     function array_compare1D(A,B) result(check)
!       real(dp), intent(in) :: A(:),B(:)
!       logical check
!       integer :: k
!       check=.false.
!       if (size(A) == size(B)) then
!          check=.true.
!          do k=1,size(A)
!             if(A(k) /= B(k) ) check=.false.
!          enddo
!       endif
!       return
!     end function array_compare1D

!     function array_compare1D0(A,B) result(check)
!       real(dp), intent(in) :: A(:),B
!       logical check
!       integer :: k
!       check=.true.
!       do k=1,size(A)
!          if(A(k) /= B ) check=.false.
!       enddo
!       return
!     end function array_compare1D0

!     function array_compare2D(A,B) result(check)
!       real(dp), intent(in) :: A(:,:),B(:,:)
!       logical check
!       integer :: k,i
!       check=.false.
!       if (size(A,1) == size(B,1) .and. size(A,2) == size(B,2) ) then
!          check=.true.
!          do k=1,size(A(1,:))
!             do i=1,size(A(:,1))
!                if(A(i,k) /= B(i,k) ) check=.false.
!             enddo
!          enddo
!       endif
!       return
!     end function array_compare2D


!     function array_compare2D0(A,B) result(check)
!       real(dp), intent(in) :: A(:,:),B
!       logical check
!       integer :: k,i
!       check=.true.
!       do k=1,size(A(1,:))
!          do i=1,size(A(:,1))
!             if(A(i,k) /= B ) check=.false.
!          enddo
!       enddo

!       return
!     end function array_compare2D0

!     function polynomial_interpolotation(x,x_array,y_array,inter_power,error) result(y)
!       ! needs a careful check...
!       ! NR 3rd edition p 118
!       ! Neville's algorithm
!       implicit none
!       real(dp), intent(in) :: x,x_array(:),y_array(:)
!       real(dp), intent(out), optional :: error
!       integer,  intent(in), optional :: inter_power
!       real(dp)  :: y
!       real(dp) :: dy,diff
!       real(dp), dimension(:), allocatable :: c,d,xi,yi
!       integer :: n,ns,i,j
!       if( size(x_array) /= size(y_array) ) stop "Incorrect size of x,y arrays in polynomial_interpolation"

!       ! select the closes point
!       diff=abs(x-x_array(1))
!       ns=1
!       do i=2,size(x_array)
!          if ( diff > abs(x-x_array(i)) ) then
!             diff = abs(x-x_array(i))
!             ns = i
!          endif
!       enddo
      
!       if (present(inter_power) ) then
!          n=min(inter_power,size(x_array))
!       else
!          n=size(x_array)
!       endif

!       allocate(xi(n),yi(n),c(n),d(n))

!       !print *,"n=",n
!       range_define: block
!         integer :: i_shift,i_start,i_stop,ns_tmp
!         i_shift=Heaviside_function(x-x_array(ns))
!         i_start=ns-n/2 + (1-2*(n/2.-n/2))*((i_shift+1)/2)
!         i_stop=ns+n/2 + (1-2*(n/2.-n/2))*((i_shift-2)/2)
!         ns_tmp= (n+1)/2 + (1-2*(n/2.-n/2))*((2-i_shift)/2.)
        
!         if ( i_start < 1 ) then
!            i_start = 1
!            i_stop = n
!            ns_tmp=ns
!         endif
!         if (i_stop > size(x_array)) then
!            i_start = size(x_array)-n+1
!            i_stop = size(x_array)
!            ns_tmp=ns+n-size(x_array)
!         endif
!         xi=x_array(i_start:i_stop)
!         yi=y_array(i_start:i_stop)
!         ns=ns_tmp
!         !        print *,"star at ",i_start,"stop at ",i_stop,"closer point ",ns,"array of x: ",(x-xi)
!         ! seems correct... oct 2017
!       end block range_define


!       ! neville: block ! this is consisten with my note, but probably there is some typo
!       !   c=yi
!       !   do i=2,n
!       !      do j=1,n-i+1
!       !         d(j)=((x-xi(j+i-1))*c(j)+(xi(j)-x)*c(j+1))/(xi(j)-xi(j+i-1))
!       !      enddo
!       !      c=d
!       !   enddo
!       !   y=c(1)
!       ! end block neville

!       nevile: block ! this is block that tracks error, consisten with my note
!         real(dp) :: ho, hp, den, w
!         integer :: k,l
!         y=yi(ns)
!         c=yi
!         d=yi
!         do k=1,n-1
!            do l=1,n-k
!               w=(c(l+1)-d(l))
!               den=xi(k+l)-xi(l)
!               ho=(x-xi(l))
!               hp=(x-xi(k+l))
!               if ( den == 0. ) then
!                  stop "error in nevile2"
!               endif
!               den= w / den
!               c(l)=ho*den
!               d(l)=hp*den
!            enddo
!            if ( 2*ns <= n-l ) then
!               dy=c(ns)
!            else
!               ns=ns-1
!               dy=d(ns)
!            endif
!            y=y+dy
!         enddo
!         if (present(error) ) then
!            if ( n /= 1) then
!               error=abs(dy)
!            else
!               error=abs(y)
!            endif
!         endif
!       end block nevile
      
!     end function polynomial_interpolotation
    


!     function array_interpolation(x,x_array,y_array,debug) result(y) ! quadratic interpolation
!       implicit none
!       real(dp), intent(in) :: x,x_array(:),y_array(:)
!       real(dp) :: y
!       logical, optional :: debug

!       integer n1,n2,n3 ! for quadratic spline

!       if( size(x_array) /= size(y_array) .or. size(x_array) < 2 ) stop "Incorrect size of x,y arrays in array_interpolation"
!       if ( x < x_array(1) .or. x>x_array(size(x_array)) ) stop "x is out of range in array_interpolation"
!       select case (size(x_array))
!       case (2)
!          if ( x_array(1) == x_array(2) )  stop "Incorrect x-array in array_interpolation"
!          y=y_array(1)+(y_array(2)-y_array(1))*(x-x_array(1))/(x_array(2)-x_array(1))
!       case default
!          n1=1
!          do while (x >= x_array(n1+1) )
!             n1=n1+1
!          enddo
!          if (n1 > 1 ) then
!             n2=n1+1
!             n3=n1-1
!          else
!             n2=n1+1 ! why these points exist???? XXXXX
!             n3=n1+2 
!          endif
!          ! this can be improved NR p 118
!          if (present(debug) ) then
!             print *,"ineterpolation", x,x_array(n1),x_array(n2),x_array(n3),y_array(n1),y_array(n2),y_array(n3)
!          endif
!          y=y_array(n1)*(x-x_array(n2))*(x-x_array(n3))/(x_array(n1)-x_array(n2))/(x_array(n1)-x_array(n3)) + &
!               &y_array(n2)*(x-x_array(n3))*(x-x_array(n1))/(x_array(n2)-x_array(n1))/(x_array(n2)-x_array(n3)) + &
!               &y_array(n3)*(x-x_array(n1))*(x-x_array(n2))/(x_array(n3)-x_array(n2))/(x_array(n3)-x_array(n1))
!       end select
!       return
!     end function array_interpolation


!     function array_interpolation_positive(x,x_array,y_array,debug) result(y) ! quadratic interpolation
!       implicit none
!       real(dp), intent(in) :: x,x_array(:),y_array(:)
!       real(dp)             :: ln_x,ln_x_array(size(x_array)),ln_y_array(size(x_array))
!       real(dp) :: y,ln_y
!       logical, optional :: debug

!       integer n1,n2 

!       if( size(x_array) /= size(y_array) .or. size(x_array) < 2 ) stop "Incorrect size of x,y arrays in array_interpolation_pos"
!       if ( x < x_array(1) .or. x > x_array(size(x_array)) ) stop "x is out of range in array_interpolation_pos"
!       ln_x = log(x)
!       ln_x_array = log(x_array) 
!       ln_y_array = log(y_array) 
!       n1 = nearest_element_scalar(ln_x_array, ln_x, -1, .true.)
!       n2 = n1 + 1
!       if ( x_array(n1) == x_array(n2) )  stop "Incorrect x-array in array_interpolation_pos"
!       ln_y = ln_y_array(n1) + &
!            & ( ln_y_array(n2) - ln_y_array(n1) ) * &
!            & ( ln_x - ln_x_array(n1) ) / (ln_x_array(n2)-ln_x_array(n1))
!       y = exp(ln_y)
!       return
!     end function array_interpolation_positive
        
    

    
!     function array_linear_interpolation1D(x,x_array,y_array) result(y)
!       implicit none
!       real(dp), intent(in) :: x,x_array(:),y_array(:)
!       real(dp) :: y
!       integer :: i
!       if( size(x_array) /= size(y_array) .or. size(x_array) < 2 ) stop "Incorrect size of x,y arrays in array_interpolation"
!       if ( x < x_array(1) .or. x>x_array(size(x_array)) ) stop "x is out of range in array_interpolation"
!       i=value_index(x_array,x)
!       if ( i == 1 ) i=2
!       y=y_array(i-1)+(x-x_array(i-1))/(x_array(i)-x_array(i-1))*(y_array(i)-y_array(i-1))
!       return
!     end function array_linear_interpolation1D
    

!     function array_linear_interpolation2D(x,y,x_array,y_array,array) result(xy_value)
!       implicit none
!       real(dp), intent(in) :: x,y,x_array(:),y_array(:),&
!            &array(size(x_array),size(y_array))
!       real(dp) :: xy_value
!       integer :: i,j
!       real(dp) :: u,t

!       i=value_index(x_array,x)
!       j=value_index(y_array,y)
!       if ( i == 1 ) i=2
!       if ( j == 1 ) j=2
!       u=(x-x_array(i-1))/(x_array(i)-x_array(i-1))
!       t=(y-y_array(j-1))/(y_array(j)-y_array(j-1))
!       xy_value=&
!            &(1._dp-u)*(1._dp-t)*array(i-1,j-1)+&
!            &u*(1._dp-t)*array(i,j-1)+&
!            &(1._dp-u)*t*array(i-1,j)+&
!            &u*t*array(i,j)
!       return
!     end function array_linear_interpolation2D


!     ! function array_linear_interpolation2D(x,x_array,y_array) result(y)
!     !   implicit none
!     !   real(dp), intent(in) :: x(2),x_array(:,:),y_array(:,:)
!     !   real(dp) :: y
!     !   integer :: i
!     !   if(    size(x_array,1) /= size(y_array,1) .or.&
!     !        & size(x_array,2) /= size(y_array,2) .or.&
!     !        & size(x_array,1) < 2 .or. &
!     !        & size(x_array,2) < 2) &
!     !        & stop "Incorrect size of x,y arrays in array_interpolation"
!     !   if ( x(1) < x_array(1,1) .or. x>x_array(size(x_array)) ) stop "x is out of range in array_interpolation"
!     !   i=value_index(x_array,x)
!     !   y=y_array(i)+(x-x_array(i))/(x_array(i+1)-x_array(i))*(y_array(i+1)-y_array(i))
!     !   return
!     ! end function array_linear_interpolation2D
    




!     subroutine spline_coeficients_v2(x,x_array,y_array,a,b,c,middle_array_index) ! diff to other: using real spline ! very time inefficient: consider using save atribute
!       ! doesn't solve the problem as compared to the fake spline (give IDENTICAL spectra)

!       ! x_array -- array of variable values
!       ! y_array -- array of function values
!       ! x -- variable value
!       ! middle_array_index -- the previous to x node of array
!       ! a,b,c quadratic spline coefficients

!       implicit none
!       real(dp), intent(in) :: x,x_array(:),y_array(:)
!       real(dp), intent(out) :: a,b,c !constant, linear, and quadratic spline coefficints
!       integer, intent(out) :: middle_array_index !location in the array
!       real(dp) :: spline(3,size(x_array)-1)
!       integer n1

!       if( size(x_array) /= size(y_array) .or. size(x_array) < 2 ) stop "Incorrect size of x,y arrays in splice_coeficients v2"
!       if ( x < x_array(1) .or. x>x_array(size(x_array)) ) stop "x is out of range in splice_coeficients v2"
!       select case (size(x_array))
!       case (2)
!          if ( x_array(1) == x_array(2) )  stop "Incorrect x-array in splice_coeficients v2"
!          a=(y_array(1)*x_array(2)-y_array(2)*x_array(1))/(x_array(2)-x_array(1))
!          b=(y_array(2)-y_array(1))/(x_array(2)-x_array(1))
!          c=0._dp
!          middle_array_index=1
!       case default
!          spline=spline_coeficient_array(x_array,y_array)
!          n1=1
!          do while (x >= x_array(n1+1) .and. n1 < size(x_array)-2  )
!             n1=n1+1
!          enddo
!          a=spline(1,n1)
!          b=spline(2,n1)
!          c=spline(3,n1)
!       end select
!       middle_array_index=n1
!       return
!     end subroutine spline_coeficients_v2


!     function spline_coeficient_array(x_array,y_array) result(spline)
!       implicit none
!       real(dp), intent(in) :: x_array(:),y_array(:)
!       real(dp) :: spline(3,size(x_array)-1)
!       integer :: i
!       real(dp) :: tmp1,tmp2

!       if (size(x_array) < 3 .or. size(x_array) /=size(y_array) ) &
!            &stop "size of arrays in not sufficient for spline extrapolation spline_coeficient_array"

!       spline(1,1)=&
!            &y_array(1)*x_array(2)*x_array(3)/(x_array(1)-x_array(2))/(x_array(1)-x_array(3))+&
!            &y_array(2)*x_array(1)*x_array(3)/(x_array(2)-x_array(1))/(x_array(2)-x_array(3))+&
!            &y_array(3)*x_array(2)*x_array(1)/(x_array(3)-x_array(2))/(x_array(3)-x_array(1))
!       spline(2,1)=-1._dp*(&
!            &y_array(1)*(x_array(2)+x_array(3))/(x_array(1)-x_array(2))/(x_array(1)-x_array(3))+&
!            &y_array(2)*(x_array(1)+x_array(3))/(x_array(2)-x_array(1))/(x_array(2)-x_array(3))+&
!            &y_array(3)*(x_array(2)+x_array(1))/(x_array(3)-x_array(2))/(x_array(3)-x_array(1)))

!       spline(3,1)=&
!            &y_array(1)/(x_array(1)-x_array(2))/(x_array(1)-x_array(3))+&
!            &y_array(2)/(x_array(2)-x_array(1))/(x_array(2)-x_array(3))+&
!            &y_array(3)/(x_array(3)-x_array(2))/(x_array(3)-x_array(1))

!       spline(:,2)=spline(:,1)
!       if ( size(x_array) > 3 ) then
!          do i=3,size(x_array)-1
!             tmp1=y_array(i)!spline(1,i-1)+x_array(i)*spline(2,i-1)+x_array(i)**2*spline(3,i-1) ! function value from left
!             tmp2=spline(2,i-1)+2._dp*x_array(i)*spline(3,i-1) ! derivative value from left
!             spline(3,i)=(y_array(i+1)-tmp1-tmp2*(x_array(i+1)-x_array(i)))/(x_array(i+1)-x_array(i))**2
!             spline(2,i)=tmp2-spline(3,i)*2._dp*x_array(i)
!             spline(1,i)=tmp1-tmp2*x_array(i)+spline(3,i)*x_array(i)**2
!          enddo
!       endif
!       return
!     end function spline_coeficient_array




!     recursive function nearest_element_recursive(A,l,u,aa,boundary) result(x)
!       ! A is a float monotonically increasing array
!       ! aa is a float value
!       ! l is lower index
!       ! u is upper index
!       ! x is the index of the element closest to aa
!       implicit none
!       real(dp), intent(in) :: A(:), aa
!       integer, optional, intent(in) :: boundary
!       integer,  intent(in) :: l, u
!       integer :: x
!       integer :: mid,b
!       if ( present(boundary) ) then
!          b = boundary
!       else
!          b = 0
!       end if
!       if ( u > l + 1 ) then
!          mid = l + ( u - l ) / 2
!          if ( A(mid) == aa ) then
!             x = mid
!          elseif ( A(mid) > aa ) then
!             x = nearest_element_recursive(A,l,mid,aa,b)
!          else
!             x = nearest_element_recursive(A,mid,u,aa,b)
!          end if
!       else
!          if ( b < 0 ) then
!             x = l
!          else if ( b > 0 ) then
!             x = u
!          else
!             if ( 2 * aa > A(u) + A(l) ) then
!                x=u
!             else
!                x=l
!             end if
!          end if
!       endif
!     end function nearest_element_recursive
    
!   function nearest_element_scalar(A, aa, boundary, strict) result(x)
!     ! A is a float monotonically increasing array
!     ! aa is a float value
!     ! boundary is an optional ingeger variable, if < 0 a lower value  if > 0, a higher value, if ==0 closest
!     ! strict is an optional logical variable, if strict  returns -1 if outside of A
!     ! x is the index of the element closest to aa
!       implicit none
!       real(dp), intent(in)     :: A(:), aa
!       integer :: x,N
!       logical, optional, intent(in) :: strict
!       integer, optional, intent(in) :: boundary
!       logical :: s
!       integer :: b
!       if ( present(boundary) ) then
!          b = boundary
!       else
!          b = 0
!       end if
!       if ( present(strict) ) then
!          s = strict
!       else
!          s = .true.
!       end if
!       N=size(A)
!       if ( aa <  A(1) .and. s ) then
!          x = -1
!       elseif  ( aa <= A(1) ) then
!          x = 1
!       elseif ( aa > A(N) .and. s ) then
!          x = -1
!       elseif ( aa >= A(N) ) then
!          x = N
!       else 
!          x = nearest_element_recursive(A, 0, N, aa, b)
!       end if
!     end function nearest_element_scalar
!   function nearest_element_vector(A, aa, boundary, strict) result(x)
!     ! A is a float monotonically increasing array
!     ! aa is a float monotonically increasing array of values
!     ! boundary is an optional ingeger variable, if < 0 a lower value  if > 0, a higher value, if ==0 closest
!     ! strict is an optional logical variable, if strict  returns -1 if outside of A
!     ! x is the index of the element closest to aa
!     implicit none
!     real(dp), intent(in)         :: A(:), aa(:)
!     integer, dimension(size(aa)) :: x
!     integer :: N, NN, i, l, b
!     logical, optional, intent(in) :: strict
!     integer, optional, intent(in) :: boundary
!     logical :: s
!     if ( present(boundary) ) then
!        b = boundary
!     else
!        b = 0
!     end if
!     if ( present(strict) ) then
!        s = strict
!     else
!        s = .true.
!     end if
!     N  = size(A)
!     NN = size(aa)
!     where ( aa <= A(1) )
!        x = 1
!     end where
!     where ( aa >= A(N) )
!        x = N
!     end where

!     if ( s ) then
!        where ( aa < A(1) .or. aa > A(N) )
!           x = -1
!        end where
!     endif
!     where ( aa > A(1) .and. aa < A(N) )
!        x = -2
!     end where
!     l = 1
!     do i=1,NN
!        if (x(i) == -2 ) then
!           x(i) = nearest_element_recursive(A, l, N, aa(i),b)
!           l = x(i)
!        end if
!     enddo
!   end function nearest_element_vector
  
    
   
    


!     function slave_integral1(a,b,c) result(int1)
!       implicit none
!       real(dp), intent(in) :: a,b,c
!       real(dp) :: int1
!       if ( a > b .or. c > a .or. c > b ) stop "slave_integral1: wrong parameters"
!       int1=sqrt(b ** 2 - c ** 2 ) - sqrt( a ** 2 - c ** 2)
!       return
!     end function slave_integral1
!     function slave_integral2(a,b,c) result(int1)
!       implicit none
!       real(dp), intent(in) :: a,b,c
!       real(dp) :: int1
!       if ( a > b .or. c > a .or. c > b ) stop "slave_integral1: wrong parameters"
!       int1=0.5 * (b * sqrt(b ** 2 - c ** 2 ) - a * sqrt( a ** 2 - c ** 2))
!       int1=int1+0.5 * c ** 2 * log( (sqrt(b **2 - c ** 2) + b )/(sqrt(a ** 2 - c ** 2) + a ))
!       return
!     end function slave_integral2
    
!     function projection_3D_to_2D(r,f,rho) result(f_int)
!       implicit none
!       real(dp), intent(in) :: r(:),f(:),rho
!       real(dp) :: f_int
!       integer :: i,k
!     ! """
!     ! r=array
!     ! f=array of values for distances r
!     ! rho=projected distance
!     ! linear interpolation for f is used
!     ! """
!       if (size(r) /= size(f)) stop "projection_3D_to_2D: array mismatch"
!       i=1
!       f_int=0.
!       if (rho < r(size(r))) then
!          do k=1,size(r)-1
!             if (rho>=r(k)) i=k+1
!          enddo
!          if ( rho < r(1) ) then
!             f_int=slave_integral1(rho,r(1),rho) * f(1)
!          else
!             f_int=slave_integral1(rho,r(i),rho) * (f(i-1)*r(i) - f(i)*r(i-1))/(r(i)-r(i-1))
!             f_int=f_int+slave_integral2(rho,r(i),rho) * (f(i) - f(i-1))/(r(i)-r(i-1))
!          endif
!          do k=i,size(r)-1
!             f_int=f_int+slave_integral1(r(k),r(k+1),rho) * (f(k)*r(k+1) - f(k+1)*r(k))/(r(k+1)-r(k))
!             f_int=f_int+slave_integral2(r(k),r(k+1),rho) * (f(k+1) - f(k))/(r(k+1)-r(k))
!          enddo
!          f_int=2.*f_int
!       endif
!     end function projection_3D_to_2D

!     subroutine trapzd(func,a,b,s,n)
!       ! NR 3rd edition page 163
!       ! trapezoidal rule
!       implicit none
!       real(dp), intent(in) :: a,b
!       integer, intent(in out) :: n
!       real(dp), intent(in out) :: s
!       interface
!          function func(x) &
!               &result(y)
!            real(8), intent(in) :: x
!            real(8) :: y
!          end function func
!       end interface
!       real(dp) :: xx,dx,sum
!       integer :: i,it_num
!       if ( n < 0 ) then
!          stop "trapzd: n<1"
!       elseif ( n == 0 ) then
!          s=0.5*(func(a)+func(b))*(b-a)
!          n=1
!       else
!          it_num=ishft(1,n-1)
!          dx=(b-a)/it_num
!          xx=a+dx/2._dp
!          sum=0._dp
!          do i=1,it_num
!             sum=sum+func(xx)
!             xx=xx+dx
!          enddo
!          n=n+1
!          s=0.5*(s+dx*sum)
!       endif
!     end subroutine trapzd
    
!     subroutine qsimp(func,a,b,relative_error,integral,result_code)
!       ! NR 3rd edition page 165
!       ! simpson integration
!       implicit none
!       real(dp), intent(in) :: a,b,relative_error
!       real(dp), intent(out) :: integral
!       integer, intent(out) :: result_code
!       interface
!          function func(x) &
!               &result(y)
!            real(8), intent(in) :: x
!            real(8) :: y
!          end function func
!       end interface
!       integer, parameter :: max_itteration=30
!       real(dp) :: s0,s,r0,r
!       integer :: n
!       result_code=0
!       n=0
!       s=0
!       s0=0
!       converge: block
!         integer :: i
!         do i=1,max_itteration
!            call trapzd(func,a,b,s,n)
!            if (i>5) then
!               r=(4.*s-s0)/3.
!               if ( abs((r-r0)/r) < relative_error .or. ( r == 0 .and. r0 == 0 ) ) then
!                  integral=r
!                  exit converge
!               endif
!               r0=r
!               s0=s
!            elseif (i==5) then
!               r0=(4.*s-s0)/3.
!               s0=s
!            elseif (i==4) then
!               s0=s
!            endif
!         end do
!         print *,"qsimp, not converged"
!         integral=(r-r)/(r-r)
!         result_code=1
!       end block converge
!     end subroutine qsimp


!     subroutine qromb(func,a,b,relative_error,integral,result_code)
!       ! NR 3rd edition page 166
!       ! Romberg integration
!       implicit none
!       real(dp), intent(in) :: a,b,relative_error
!       real(dp), intent(out) :: integral
!       integer, intent(out) :: result_code
!       interface
!          function func(x) &
!               &result(y)
!            real(8), intent(in) :: x
!            real(8) :: y
!          end function func
!       end interface
!       integer, parameter :: max_itteration=20
!       integer, parameter :: interpolation_power=5
!       real(dp), dimension(max_itteration) :: h,ff
!       real(dp) :: s,sromb,error
!       integer :: j,i,n
!       n=0
!       h(1)=1.
!       converge: block
!       do j=1,max_itteration
!          call trapzd(func,a,b,s,n)
!          ff(j)=s
!          if (j >= interpolation_power) then
!             sromb=polynomial_interpolotation(&
!                  &0._dp,h(j-interpolation_power+1:j),&
!                  &ff(j-interpolation_power+1:j),&
!                  &interpolation_power,error)
!             if ( abs(error) < abs(sromb*relative_error)) then
!                integral=sromb
!                exit converge
!             endif
            
!          endif
!          if ((j+1) <= max_itteration ) h(j+1)=0.25*h(j)
!       enddo
!       print *,"romb integration, not converged"

!     end block converge
!   end subroutine qromb
    
! !!! "System tools"
!   function to_upper_case(strIn) result(strOut)
!     ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
!     ! Original author: Clive Page

!     implicit none

!     character(len=*), intent(in) :: strIn
!     character(len=len(strIn)) :: strOut
!     integer :: i,j

!     do i = 1, len(strIn)
!        j = iachar(strIn(i:i))
!        if (j>= iachar("a") .and. j<=iachar("z") ) then
!           strOut(i:i) = achar(iachar(strIn(i:i))-32)
!        else
!           strOut(i:i) = strIn(i:i)
!        end if
!     end do
!     strOut = TRIM(strOut)
!   end function to_upper_case


! ! !!!! General purpose statistical subroutines 
!     function gauss_distribution(mu, sigma) result(z)
!       ! Box–Muller
!       implicit none
!       real(dp), intent(in) :: mu, sigma
!       real(dp)             :: z
!       real(dp)             :: x(2), s, z0
!       do
!          call random_number(x)
!          x = (1. - 2. * x)
!          s = x(1) ** 2 + x(2) ** 2
!          if ( s > 0 .and. s <= 1) then
!             z0 = x(1) * sqrt(-2 * log(s) / s)
!             z = mu + sigma * z0
!             exit
!          endif
!       enddo
!     end function gauss_distribution

! !!!! General purpose physical subroutines 

!     function stellar_radius(L_s,T_s) result(R_s) ! all in cgs units
!       implicit none
!       real(dp), intent(in) :: L_s,T_s
!       real(dp) :: R_s
!       R_s=sqrt(L_s/(4.*pi*T_s**4*stefan_boltzmann))
!       return
!     end function stellar_radius


!     function stellar_luminosity(R_s,T_s) result(L_s) ! all in cgs units
!       implicit none
!       real(dp), intent(in) :: R_s,T_s
!       real(dp) :: L_s
!       L_s=4.*pi*T_s**4*stefan_boltzmann*R_s**2
!       return
!     end function stellar_luminosity

!     ! an attempt to do something smart with obtional arguments
!     ! function stellar_parameters(radius,luminosity,temperature) ! all in cgs units
!     !   implicit none
!     !   real(dp), optional :: radius,luminosity,temperature
!     !   real(dp) :: tmp
!     !   if ( .not. present(radius) ) then
!     !      if (present(luminosity) .and. present(temperature)) then
!     !         tmp=sqrt(luminosity/(4.*pi*temperature**4*stefan_boltzmann))
!     !      else
!     !         stop "not enough parameters in stellar_parameters"
!     !      endif
!     !   endif
!     !   if ( .not. present(luminosity) ) then
!     !      if (present(radius) .and. present(temperature)) then
!     !         tmp=4.*pi*temperature**4*stefan_boltzmann*radius**2
!     !      else
!     !         stop "not enough parameters in stellar_parameters"
!     !      endif
!     !   endif
!     ! end function stellar_parameters
    
    

!     ! two reference frames, a photon moves in the rest system at angle cos_angle_rest_frame to
!     ! the velocity of the moving system. In the moving systems it has angle cos_angle_moving_frame
!     ! to the direction oposide to the systems' relative velocity (in the moving system!)
!     elemental function photon_abberation(cos_angle_rest_frame,frame_velocity) &
!          &result(cos_angle_moving_frame)
!       implicit none
!       real(dp), intent(in) :: cos_angle_rest_frame,frame_velocity
!       real(dp) :: cos_angle_moving_frame
!       real(dp) :: tmp
!       tmp=(1._dp-frame_velocity*cos_angle_rest_frame)
!       cos_angle_moving_frame=(cos_angle_rest_frame-frame_velocity)/tmp
!       return
!     end function photon_abberation




!     ! two reference frames, a photon moves in the rest system at angle cos_angle_rest_frame to
!     ! the velocity of the moving system. In the moving systems it has energy energy_moving_frame
!     elemental function photon_energy_relativistic_transformation&
!          &(energy_rest_frame,cos_angle_rest_frame,frame_velocity) result(energy_moving_frame)
!       implicit none
!       real(dp), intent(in) :: energy_rest_frame,cos_angle_rest_frame,frame_velocity
!       real(dp) :: energy_moving_frame
!       real(dp) :: frame_relative_lorentz_factor,tmp
!       frame_relative_lorentz_factor=1._dp/sqrt(1._dp-frame_velocity**2)
!       energy_moving_frame=energy_rest_frame*&
!            &(1._dp-frame_velocity*cos_angle_rest_frame)*frame_relative_lorentz_factor
!       return
!     end function photon_energy_relativistic_transformation

    
!     !   Photons density (photon/cm^3) for Plank distribution. Ephoton in mc^2,
!     !   Temperature Tmcc in mc^2
!     subroutine Planck(photon_energy,temperature_mcc,photon_density)
!       implicit none
!       real(dp), intent(in) ::photon_energy,temperature_mcc
!       real(dp), intent(out) ::photon_density
!       photon_density=Planck_function(photon_energy,temperature_mcc)
!       return
!     end subroutine Planck


!     elemental function Planck_function(photon_energy,temperature_mcc) result(photon_density)
!       implicit none
!       real(dp), intent(in) ::photon_energy,temperature_mcc
!       real(dp) photon_density

!       real(dp) :: t3,t2,t6,t7
!       if (temperature_mcc > 0._dp) then
!          t3 = photon_energy**2
!          t2=photon_energy/temperature_mcc
!          t6=exp(-t2)
!          t7=1._dp-t6
!          if( t2 < 1.e-4_dp) t7=t2
!          photon_density = 0.17595e31_dp*t3*t6/t7 ! m_e^3 c^3 / (pi^2 \hbar^3)
!       else
!          photon_density=0._dp
!       endif

!       return
!     end function Planck_function

    

!     subroutine greybody(photon_energy,temperature_mcc,field_energy_density_mcc_per_cm3,photon_density)
!       implicit none
!       real(dp), intent(in) ::photon_energy,temperature_mcc,field_energy_density_mcc_per_cm3
!       real(dp), intent(out) :: photon_density
!       photon_density=greybody_function(photon_energy,temperature_mcc,field_energy_density_mcc_per_cm3)
!       return
!     end subroutine greybody
    
!     elemental function greybody_function(photon_energy,temperature_mcc,field_energy_density_mcc_per_cm3) result(photon_density)
!       implicit none
!       real(dp), intent(in) ::photon_energy,temperature_mcc,field_energy_density_mcc_per_cm3
!       real(dp) photon_density

!       real(dp) :: t1
!       if (field_energy_density_mcc_per_cm3 > 0._dp) then
!          t1 = field_energy_density_mcc_per_cm3/(1.142643e+31_dp*temperature_mcc**4)
!          photon_density = t1*Planck_function(photon_energy,temperature_mcc)
!       else
!          photon_density=0._dp
!       endif

!       return
!     end function greybody_function
    
!     function greybody_dilution(temperature_mcc,field_energy_density_mcc_per_cm3) result(dilution_coef)
!       implicit none
!       real(dp), intent(in) ::temperature_mcc,field_energy_density_mcc_per_cm3
!       real(dp) dilution_coef

!       if (field_energy_density_mcc_per_cm3 > 0._dp .and. temperature_mcc > 0._dp ) then
!          dilution_coef = field_energy_density_mcc_per_cm3/(1.142643e+31_dp*temperature_mcc**4)
!       else
!          dilution_coef=0._dp
!       endif

!       return
!     end function greybody_dilution
    

!     ! Keplerian motion
!     function periastron_separation(orbit_eccentricity,major_semiaxis) result(periastron)
!       implicit none 
!       real(dp), intent(in) :: orbit_eccentricity,major_semiaxis
!       real(dp) :: periastron
!       periastron=major_semiaxis*(1.-orbit_eccentricity)
!     end function periastron_separation

!     function apastron_separation(orbit_eccentricity,major_semiaxis) result(apastron)
!       implicit none 
!       real(dp), intent(in) :: orbit_eccentricity,major_semiaxis
!       real(dp) :: apastron
!       apastron=major_semiaxis*(1.+orbit_eccentricity)
!     end function apastron_separation
    

!     function major_semiaxis(period_days,system_mass) result(semiaxis_cm)
!       implicit none
!       real(dp), intent(in) :: period_days,system_mass
!       real(dp) :: semiaxis_cm
      
!       semiaxis_cm=(period_days*s_per_day/(2*pi)*sqrt(gravitational_contant*system_mass*sun_mass))**0.6666
!     end function major_semiaxis
    
!     ! COORDINATE SYSTEMS
!     function to_decart(r0,optional_cor_system) result(r_dec)
!       implicit none
!       real(dp), intent(in), dimension(3) :: r0
!       real(dp), dimension(3) :: r_dec
!       integer, intent(in), optional :: optional_cor_system
!       integer :: cor_system
!       cor_system = 1
!       if ( present(optional_cor_system) ) then
!          cor_system = optional_cor_system
!       end if
!       select case (cor_system)
!       case (1)
!          r_dec(1) = r0(1) * sin(r0(2)) * cos(r0(3))
!          r_dec(2) = r0(1) * sin(r0(2)) * sin(r0(3))
!          r_dec(3) = r0(1) * cos(r0(2))
!       case default
!          r_dec=r0
!          print *,"wrong coordinate system in to_decart"
!       end select
      
!     end function to_decart

!     function to_spherical(r0,optional_cor_system) result(r_sph)
!       ! (r,theta,phi)
!       implicit none
!       real(dp), intent(in), dimension(3) :: r0
!       real(dp), dimension(3) :: r_sph
!       integer, intent(in), optional :: optional_cor_system
!       integer :: cor_system
!       cor_system = 1
!       if ( present(optional_cor_system) ) then
!          cor_system = optional_cor_system
!       end if
!       select case (cor_system)
!       case (1)
!          r_sph(1) = sqrt(sum(r0**2))
!          r_sph(2) = 0.
!          if ( r_sph(1) /= 0 ) then
!             r_sph(2) = acos( r0(3) / r_sph(1) )
!          endif
!          r_sph(3) = atan2( r0(2), r0(1) ) 
!       case default
!          r_sph=r0
!          print *,"wrong coordinate system in to_spherical"
!       end select
      
!     end function to_spherical
    
    
!     function euler_rotation(phi,theta,psi,r1) result(r2)
!       ! transforms to decart coordinated system with orientation determined by the Euler angles
!       implicit none
!       real(dp), intent(in) :: phi,theta,psi ! notation as in LL1par35
!       real(dp), intent(in), dimension(3) :: r1 ! coordinate  of a vector: i.e., system (XYZ) as defined in LL1par35
!       real(dp), dimension(3) :: r2 ! coordinated of the vector in the rotated system: i.e, (x1,x2,x3) as defined in LL1par35

!       ! in for the standard orbital elements in binary systems Z is line sight (i.e. XY is the plane of sky)
!       !                                                        x3 is selected by the rotation
!       !                                                        X is node line (ON as defined in LL1par35), i.e., phi=0
!       !                                                        theta is orbital inclination i
!       !                                                        psi is longitude of periastron (x1 is directed to periastron)

!       real(dp), dimension(3,3) :: euler_matrix

!       ! if (  phi<0. .or. phi>2.*pi .or. &
!       !      &psi < 0. .or. psi > 2*pi .or. &
!       !      &theta<0. .or. theta > pi ) stop "non-standard Euler angles values"
!       if ( theta<0. .or. theta > pi ) stop "non-standard Euler angles values"

!       ! Fortran fils the matrixes by column (from the most left index)
!       euler_matrix=reshape(&
!            &[cos(phi)*cos(psi)-sin(phi)*sin(psi)*cos(theta),-sin(phi)*cos(psi)*cos(theta)-cos(phi)*sin(psi),sin(phi)*sin(theta),&
!            &cos(phi)*sin(psi)*cos(theta)+sin(phi)*cos(psi),cos(phi)*cos(psi)*cos(theta)-sin(phi)*sin(psi),-cos(phi)*sin(theta),&
!            &sin(psi)*sin(theta),cos(psi)*sin(theta),cos(theta)],[3,3])


!       r2=matmul(euler_matrix,r1)

!     end function euler_rotation

!     function reverse_euler_rotation(phi,theta,psi,r1) result(r2)
!       implicit none
!       real(dp), intent(in) :: phi,theta,psi ! notation as in LL1par35
!       real(dp), intent(in), dimension(3) :: r1 ! coordinate  of a vector in the rotated system: i.e, (x1,x2,x3) as defined in LL1par35
!       real(dp), dimension(3) :: r2 ! coordinated of the vector: i.e., system (XYZ) as defined in LL1par35

!       ! in for the standard orbital elements in binary systems Z is line sight (i.e. XY is the plane of sky)
!       !                                                        x3 is selected by the rotation
!       !                                                        X is node line (ON as defined in LL1par35), i.e., phi=0
!       !                                                        theta is orbital inclination i
!       !                                                        psi is longitude of periastron (x1 is directed to periastron)

!       real(dp), dimension(3,3) :: euler_matrix

!       ! if (  phi<0. .or. phi>2.*pi .or. &
!       !      &psi < 0. .or. psi > 2*pi .or. &
!       !      &theta<0. .or. theta > pi ) stop "non-standard Euler angles values"
!       if ( theta<0. .or. theta > pi ) stop "non-standard Euler angles values"

!       ! Fortran fils the matrixes by column (from the most left index) reverse matrix
!       euler_matrix=reshape(&
!            &[cos(phi)*cos(psi)-sin(phi)*sin(psi)*cos(theta),cos(phi)*sin(psi)*cos(theta)+sin(phi)*cos(psi),sin(psi)*sin(theta),&
!            &-sin(phi)*cos(psi)*cos(theta)-cos(phi)*sin(psi),cos(phi)*cos(psi)*cos(theta)-sin(phi)*sin(psi),cos(psi)*sin(theta)&
!            &,sin(phi)*sin(theta),-cos(phi)*sin(theta),cos(theta)],[3,3])


!       r2=matmul(euler_matrix,r1)

!     end function reverse_euler_rotation
    


end module utils


  
