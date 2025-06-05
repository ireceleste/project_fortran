module FunctionModule
    use utils
    use Minimizer
    use HistogramHandler

    implicit none

contains

     function FCN(fitter_instance) result(val)
            
            type(Fitter) :: fitter_instance
            real(dp) :: val

            select case(fitter_instance%likelihood)
        
                
            case(1)  ! Unbinned Fit
                if(fitter_instance%save_pdf) then
                    val = chi2binnedFCN(fitter_instance%model, fitter_instance%yvals, fitter_instance%err_yvals, &
                                            fitter_instance%xvals, fitter_instance%p_min, &
                                            fitter_instance%pdf_vals)
                else
                    val = chi2binnedFCN(fitter_instance%model, fitter_instance%yvals, fitter_instance%err_yvals, &
                                            fitter_instance%xvals, fitter_instance%p_min)
                end if
            
            case(2,3)  ! Binned fit

                if(fitter_instance%save_pdf) then
                    val = histbinnedFCN(fitter_instance%model, fitter_instance%likelihood ,fitter_instance%hist, fitter_instance%p_min, fitter_instance%pdf_vals)
                else
                    val = histbinnedFCN(fitter_instance%model, fitter_instance%likelihood, fitter_instance%hist, fitter_instance%p_min)
                end if

            

            case default
                val = 0.0_dp  
                print *, "Warning: Likelihood index not recognized in FCN. Returning 0."
                print *, "Index:", fitter_instance%likelihood
            end select
        
    end function FCN
    
   
    function gaussian1D(x, mean, stddev) result(val)
        real(dp), intent(in) :: x, mean, stddev
        real(dp) :: val
        if(stddev <= 0.0_dp) then
            ! print *, "Error: Standard deviation must be positive."
            val = 0.0_dp 
        else 
            val = exp(-0.5_dp * ((x - mean) / stddev)**2) / (stddev * sqrt(2.0_dp * pi))
        end if
    end function gaussian1D

    function exponential1D(x, a) result(val)
        real(dp), intent(in) :: x, a
        real(dp) :: val
        if (a <= 0.0_dp) then
            ! print *, "Error: Exponential decay rate must be positive."
            val = 0.0_dp 
        else
            val = a * exp(-a * x)
        end if
    end function exponential1D

    ! Function to calculate the chi-squared value for a histogram fit to a Gaussian function

    function histbinnedFCN(model, likelihood, hist, pars, predictions) result(val)
        character(len=*), intent(in) :: model
        integer, intent(in) :: likelihood
        type(Histogram), intent(in) :: hist
        real(dp), intent(in) :: pars(:)
        real(dp) :: val, prediction
        real(dp), dimension(:), intent(inout), optional :: predictions
        integer :: i

        val = 0.0_dp

        do i = 1, hist%nbins

            if (hist%bin_errors(i)  > 0.0_dp) then  
            
                if(model == 'Gauss') then
                    prediction = pars(3) * hist%bin_widths(i) * gaussian1D(hist%bin_centers(i), pars(1), pars(2))
                else if(model == 'expo') then
                    prediction = pars(2) * hist%bin_widths(i) * exponential1D(hist%bin_centers(i), pars(1))
                else
                    print *, "Error: Invalid binned model type: ", model
                    return
                end if
                
                if (present(predictions)) then
                    predictions(i) = prediction
                end if

                if(likelihood == 2) then 
                    if(prediction > 0) val = val + 2.0_dp * (prediction - hist%bin_counts(i) + hist%bin_counts(i) * log(hist%bin_counts(i) / prediction))

                else if(likelihood == 3) then 
                    if(hist%bin_errors(i) > 0) val = val + ((hist%bin_counts(i) - prediction) / hist%bin_errors(i))**2
                else
                    print *, "Error: Invalid likelihood type in histbinnedFCN."
                    return
                end if
            end if
        end do
        
    end function histbinnedFCN




    ! Function to calculate the chi-squared value for a linear model fit to a dataset

    function chi2binnedFCN(model, yvals, err_yvals, xvals, pars, predictions) result(val)
        character(len=*), intent(in) :: model
        real(dp), intent(in) :: yvals(:), err_yvals(:), xvals(:), pars(:)
        real(dp), dimension(:), intent(inout), optional :: predictions
        real(dp) :: val, prediction
        integer :: i

        
        val = 0.0_dp

        do i = 1, size(yvals)

            if(model == 'linear') then
                prediction = pars(1) * xvals(i) + pars(2)
            else if(model == 'quadratic') then
                prediction = pars(1) * xvals(i)**2 + pars(2) * xvals(i) + pars(3)
            else
                print *, "Error: Invalid binned model type: ", model
                return
            end if

            if (present(predictions)) then
                predictions(i) = prediction
            end if

            if (err_yvals(i) > 0.0_dp) then  
                val = val + ((yvals(i) - prediction) / err_yvals(i))**2
            end if
        end do
        
    end function chi2binnedFCN

    ! Function to calculate the chi-squared value for a quadratic model fit to a dataset

    

end module FunctionModule

