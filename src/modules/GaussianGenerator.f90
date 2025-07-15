!!!!!!!!!!!!!!!!!!!!!!!!! GaussianGenerator.f08 !!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains subroutines for generating Gaussian distributions and datasets with a fixed pdf
  ! generate_gaussian: generates a Gaussian distribution using the Box-Muller transform
  ! generate_dataset_pdf: generates a dataset based on a given pdf using Monte Carlo acceptance-rejection (exponential or Gaussian)
  ! generate_dataset_model: generates a dataset based on a given model with Gaussian noise
! Called by the InputOutput module

module GaussianGenerator
    use utils
    implicit none

contains

  ! This subroutine generates a Gaussian distribution with given mean and standard deviation
  ! using the Box-Muller transform method: https://doi.org/10.1214/aoms/1177706645
  

  subroutine generate_gaussian(mu, sigma, n, output_yvals)
    real(dp), intent(in) :: mu, sigma
    integer, intent(in) :: n
    real(dp), intent(out) :: output_yvals(n)
    real(dp) :: u1, u2, r, theta, z0, z1
    integer :: i
    integer, allocatable :: seed(:)
    integer :: seed_size

    if (sigma <= 0.0_dp) then
      write(*,*) "Error: Standard deviation sigma must be positive."
      stop
    end if    
    
    call random_seed()

    i = 1
    do while (i <= n)
      call random_number(u1)
      call random_number(u2)

      if (u1 == 0.0_dp) cycle  ! Avoid log(0)
      r = sqrt(-2.0_dp * log(u1))
      theta = 2.0_dp * pi * u2

      z0 = r * cos(theta)
      z1 = r * sin(theta)

      output_yvals(i) = mu + sigma * z0
      if (i < n) then
        output_yvals(i+1) = mu + sigma * z1
        i = i + 2
      else
        i = i + 1
      end if
    end do
  end subroutine generate_gaussian


  ! This subroutine generates a dataset based on a given pdf based on the Monte Carlo
  ! acceptance-rejection sampling

  subroutine generate_dataset_pdf(n, output_yvals, pdf , pars)
    implicit none

    integer, intent(in) :: n
    character(len=*), intent(in) :: pdf
    real(dp), intent(in) :: pars(:)
    real(dp), dimension(n), intent(out) :: output_yvals

    real(dp) :: x, y, fx, xmin, xmax, fmax
    integer :: count, attempts
    real(dp) :: rng1, rng2

    integer, allocatable :: seed(:)
    integer :: seed_size

    call random_seed()

    select case (trim(pdf))

    case ("expo")
        if(size(pars) /= 1) then
            write(*, '( "Exponential PDF requires one parameter (lambda)" )')
            stop
        end if

        if(pars(1) <= 0.0_dp) then
            write(*, '( "Error: Exponential decay rate must be positive." )')
            stop
        end if

        xmin = 0.0_dp
        xmax = 10.0_dp / pars(1)
        fmax = pars(1) 

    case ("Gauss")
        if(size(pars) /= 2) then
            write(*, '( "Gaussian PDF requires two parameters (mean, sigma)." )')
            stop
        end if

        if(pars(2) <= 0.0_dp) then
            write(*, '( "Error: Standard deviation must be positive." )')
            stop
        end if
        
        xmin = pars(1) - 10.0_dp * pars(2)
        xmax = pars(1) + 10.0_dp * pars(2)
        fmax = (1.0_dp / (pars(2) * sqrt(2.0_dp * pi))) 
    case default
         write(*, '( "Unknown PDF: ", A )') trim(pdf)
        stop
    end select

    ! Acceptance-rejection sampling

    count = 0
    attempts = 0

    do while (count < n .and. attempts < 10*n)
      
      call random_number(rng1)
      call random_number(rng2)
      x = xmin + rng1 * (xmax - xmin)
      y = rng2 

      select case (trim(pdf))

      case ("expo")
          fx = pars(1) * exp(-pars(1) * x)
      case ("Gauss")
          fx = (1.0_dp / (pars(2) * sqrt(2.0_dp * pi))) * exp(-0.5_dp * ((x - pars(1)) / pars(2))**2)
      end select
      
      if (y < fx/fmax) then
          output_yvals(count+1) = x
          count = count + 1
      end if

      attempts = attempts + 1

    end do

    if (count < n) write(*, '("Warning: only generated ", I0,  " points out of ", I0)')  count, n
    

  end subroutine generate_dataset_pdf

  ! This subroutine generates a set of points (x, y, sigmay) based on a given prediction 
  ! with a Gaussian noise

  subroutine generate_dataset_model(n, output_yvals, pars, model, xvals, yerrs)
    integer, intent(in) :: n
    real(dp), intent(out) :: output_yvals(n)
    real(dp), intent(in) :: pars(:)
    character(len=*), intent(in) :: model
    real(dp), intent(out), optional :: xvals(n), yerrs(n)
    

    real(dp), dimension(n) :: random_shifts
    integer :: i
    real(dp) :: x, prediction, err

    call generate_gaussian(0.0_dp, 1.0_dp, n, random_shifts)

    do i = 1, n
      
      if (present(xvals)) then
        x = xvals(i)
      else
        x = real(i, dp)  
      end if
     
      if (present(yerrs)) then
        err = yerrs(i)
      else
        err = 1.0_dp  
      end if

      ! Model prediction to estimate the y value

      select case (model)
      case ("linear")
          if (size(pars) < 2) then
              write(*,*) "Error: linear model requires at least 2 parameters (slope, intercept)"
              stop
          end if
          prediction = pars(1) * x + pars(2)

      case ("quadratic")
          if (size(pars) < 3) then
              write(*,*) "Error: quadratic model requires 3 parameters (a, b, c)"
              stop
          end if
          prediction = pars(1) * x**2 + pars(2) * x + pars(3)

      case default
          write(*,*) "Unknown model type: ", trim(model)
          stop
      end select

      
      output_yvals(i) = prediction + random_shifts(i)*err

    end do
  end subroutine
  
end module GaussianGenerator
