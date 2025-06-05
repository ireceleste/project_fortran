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


    ! Fix the seed
    
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 12349
    call random_seed(put=seed)
    ! call random_seed()

    

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


  ! This subroutine generates a dataset based on a given pdf

  subroutine generate_dataset_pdf(n, output_yvals, pdf , pars)
    implicit none

    integer, intent(in) :: n
    character(len=*), intent(in) :: pdf
    real(dp), intent(in) :: pars(:)
    real(dp), dimension(n), intent(out) :: output_yvals

    real(dp) :: x, y, fx, xmin, xmax
    integer :: i, count
    logical :: accepted
    real(dp) :: rng1, rng2

    integer, allocatable :: seed(:)
    integer :: seed_size

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 12349
    call random_seed(put=seed)

    ! call random_seed()

    


    select case (trim(pdf))

    case ("expo")
        if(size(pars) /= 1) then
            write(*, '( "Exponential PDF requires one parameter (lambda)" )')
            stop
        end if
        xmin = 0.0_dp
        xmax = 10.0_dp * pars(1)

    case ("Gauss")
        if(size(pars) /= 2) then
            write(*, '( "Gaussian PDF requires two parameters (mean, sigma)." )')
            stop
        end if
        
        xmin = pars(1) - 10.0_dp * pars(2)
        xmax = pars(1) + 10.0_dp * pars(2)
    case default
         write(*, '( "Unknown PDF: ", A )') trim(pdf)
        stop
    end select

    
    count = 0

    do i = 1, 10 * n
        accepted = .false.

        do while (.not. accepted .and. count < n)
        
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
            
            if (y < fx) then
                output_yvals(i) = x
                count = count + 1
                accepted = .true.
            end if

        end do
    end do
end subroutine generate_dataset_pdf

  ! This subroutine generates a dataset based on a fiven model assuming a Gaussian distribution

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

      if (model == "linear") then
        prediction = pars(1) * x + pars(2)
      else if (model == "quadratic") then
        prediction = pars(1) * x**2 + pars(2) * x + pars(3)
      else if (model == "exponential") then
        prediction = pars(1) * exp(pars(2) * x)
      else
        write(*,*) "Unknown model type"
        return
      end if      
      
      output_yvals(i) = prediction + random_shifts(i)*err

    end do
  end subroutine
  
end module GaussianGenerator
