

module Minimizer
    use utils, only: dp
    use FitterModule
    use HistogramHandler
    

contains


    ! Minimize the function FCN using gradient descent

    subroutine gradient_minimize(fitter_instance, FCN, verbose_log) 
        

        interface
            function FCN(fitter_instance) result(val)
                import :: dp
                import :: Fitter
                type(Fitter) :: fitter_instance
                real(dp) :: val
            end function FCN
        end interface

        type(Fitter), intent(inout) :: fitter_instance
        character(len=*), intent(in), optional :: verbose_log
        real(dp), dimension(fitter_instance%npars) :: p, grad
        real(dp) :: h, f1, f2, fmin
        integer :: i, j, iu

        if(present(verbose_log)) then
            iu = 99
            open(unit=iu, file=verbose_log, status='replace')
            write(iu, '(/,A,/)') repeat('=', 90)
            write(iu, '("Starting gradient descent minimization with the following settings: " /)')
            write(iu, '("Likelihood: ", I0)') fitter_instance%likelihood
            write(iu, '("Model: ", A)') fitter_instance%model
            write(iu, '("Parameters: ", I0)') fitter_instance%npars
            write(iu, '("Initial parameters: ")', advance='no')
            do j = 1, fitter_instance%npars
                write(iu, '(F0.3)', advance='no') fitter_instance%p0(j)
                if (j < fitter_instance%npars) then
                    write(iu, '(A)', advance='no') ", "
                else
                    write(iu, '(/)', advance='no') 
                end if
            end do
            
            write(iu, '("Step sizes: ")', advance='no')
            do j = 1, fitter_instance%npars
                write(iu, '(ES0.2)', advance='no') fitter_instance%step(j)
                if (j < fitter_instance%npars) then
                    write(iu, '(A)', advance='no') ", "
                else
                    write(iu, '(/)', advance='no') 
                end if
            end do
            write(iu, '("Learning rates: ")', advance='no')
            do j = 1, fitter_instance%npars
                write(iu, '(ES0.2)', advance='no') fitter_instance%lrate(j)
                if (j < fitter_instance%npars) then
                    write(iu, '(A)', advance='no') ", "
                else
                    write(iu, '(/)', advance='no') 
                end if
            end do
            write(iu, '("Tolerance: ", ES0.2)') fitter_instance%tol
            write(iu, '("Maximum iterations: ", I0)') fitter_instance%max_iter
            write(iu, '(/,A,/)') repeat('=', 90)
        end if

        p = fitter_instance%p0
        h = 1.0e-6_dp  ! Step size for numerical gradient

        do i = 1, fitter_instance%max_iter 

            grad = 0.0_dp
            fmin = FCN(fitter_instance) ! Evaluate the function at the current parameters

            do j = 1, fitter_instance%npars
                
                ! Compute partial derivative with respect to p(j)
                if(fitter_instance%step(j) == 0.0_dp) then
                    print *, "Error: Step size is zero for parameter ", j
                    grad(j) = 0.0_dp
                else
                    fitter_instance%p_min(j) = p(j) + fitter_instance%step(j)
                    f1 = FCN(fitter_instance)
                    fitter_instance%p_min(j) = p(j) - fitter_instance%step(j)
                    f2 = FCN(fitter_instance)
                    grad(j) = (f1 - f2) / (2*fitter_instance%step(j))

                    fitter_instance%p_min(j) = p(j)  ! Reset to original value
                end if

                if(present(verbose_log)) then
                    write(iu, '("Iteration " , I0, " - parameter " , I0 , " = " , F11.3, " - FCN = ", F11.3, " - grad = ", ES12.3)') i,j, p(j), fmin, grad(j)
                end if
            end do
            
            if (maxval(abs(grad)) < fitter_instance%tol) exit
            p = p - fitter_instance%lrate * grad
            fitter_instance%p_min = p  
            
        end do

        fitter_instance%niter = i
        fitter_instance%p_min = p
        fitter_instance%grad_p_min = grad
        fitter_instance%min_FCN = fmin
        
        if(present(verbose_log)) then
            write(iu, '("Finished gradient descent minimization after ", I0, " / " , I0 , " iterations.")') fitter_instance%niter, fitter_instance%max_iter
            write(iu, '("Minimum function value = ", F20.5)') fmin
            write(iu, '("Minimum parameters: ")', advance='no')
            do j = 1, fitter_instance%npars
                write(iu, '(F12.3)', advance='no') fitter_instance%p_min(j)
                if (j < fitter_instance%npars) write(iu, '(A,/)', advance='no') ", "
            end do
            write(iu, '("Gradient at minimum: ")', advance='no')
            do j = 1, fitter_instance%npars
                write(iu, '(F22.3)', advance='no') fitter_instance%grad_p_min(j)
                if (j < fitter_instance%npars) write(iu, '(A,/)', advance='no') ", "
            end do
            write(iu, '("Function value at minimum = ", F20.5 , " / ndof ", I0)') fmin, fitter_instance%ndof
            write(iu, '( 1x  )' )
            close(iu)
        end if


        
    end subroutine gradient_minimize

    
    ! Estimate the covariance matrix using the Hessian

    subroutine estimate_covariance(fitter_instance, FCN)
        
        interface
            function FCN(fitter_instance) result(val)
                import :: dp
                import :: Fitter
                type(Fitter) :: fitter_instance
                real(dp) :: val
            end function FCN
        end interface

        type(Fitter), intent(inout) :: fitter_instance

        real(dp), allocatable :: cov_matrix(:,:)

        integer :: i, j, n
        real(dp) :: f0, fpp, fpm, fmp, fmm, h
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: Hess(:,:)
        type(Fitter) :: tmp

        n = fitter_instance%npars
        allocate(Hess(n, n), cov_matrix(n,n), p(n))
        p = fitter_instance%p_min

        f0 = FCN(fitter_instance)

        do i = 1, n
            do j = i, n

                if(fitter_instance%step(i) == 0.0_dp .or. fitter_instance%step(j) == 0.0_dp) then
                    Hess(i,j) = 0.0_dp
                    cycle
                end if

            
                ! fpp: +i, +j
                tmp = fitter_instance
                tmp%p_min = p
                tmp%p_min(i) = p(i) + fitter_instance%step(i)
                tmp%p_min(j) = p(j) + fitter_instance%step(j)
                fpp = FCN(tmp)

                ! fmm: -i, -j
                tmp = fitter_instance
                tmp%p_min = p
                tmp%p_min(i) = p(i) - fitter_instance%step(i)
                tmp%p_min(j) = p(j) - fitter_instance%step(j)
                fmm = FCN(tmp)



                if( i == j ) then
                    Hess(i,j) = (fpp - 2.0_dp * f0 + fmm) / (fitter_instance%step(i)**2)
                else
                    ! fpm: +i, -j
                    tmp = fitter_instance
                    tmp%p_min = p
                    tmp%p_min(i) = p(i) + fitter_instance%step(i)
                    tmp%p_min(j) = p(j) - fitter_instance%step(j)
                    fpm = FCN(tmp)

                    ! fmp: -i, +j
                    tmp = fitter_instance
                    tmp%p_min = p
                    tmp%p_min(i) = p(i) - fitter_instance%step(i)
                    tmp%p_min(j) = p(j) + fitter_instance%step(j)
                    fmp = FCN(tmp)
                    Hess(i,j) = (fpp - fpm - fmp + fmm) / (4*fitter_instance%step(i)* fitter_instance%step(j))
                end if
            end do
        end do

        ! Symmetrize Hessian
        do i = 1, n
            do j = i+1, n
                Hess(j,i) = Hess(i,j)
            end do
        end do

        
        ! Invert the Hessian to get the covariance matrix

        call invert_matrix(Hess, cov_matrix)
        fitter_instance%cov_matrix = cov_matrix

        ! Calculate the standard deviations from the covariance matrix
        do i = 1, n
            fitter_instance%p_err(i) = sqrt(cov_matrix(i,i))
        end do

    end subroutine estimate_covariance

    

    ! Print the fit results to a file

    subroutine print_fit_results(fitter_instance, filename)

        use iso_fortran_env, only: output_unit

        type(Fitter), intent(in) :: fitter_instance
        character(len=*), intent(in), optional :: filename
        integer :: iu
        character(len=8)  :: date_str
        character(len=10) :: time_str

        if(present(filename)) then
            iu = 99
            open(unit=iu, file=filename, status='replace')
        else
            iu = output_unit  
        end if


       
        write(iu, '(/,A)', advance = 'no') repeat('=', 34)
        write(iu, '(A)', advance = 'no') " FIT RESULTS SUMMARY ="
        write(iu, '(A,/)') repeat('=', 34)

        

        call date_and_time(date=date_str, time=time_str)

        write(iu, '(A,1x,A," at ",A)') "Fit performed on", date_str(1:4)//"-"//date_str(5:6)//"-"//date_str(7:8), &
                                            time_str(1:2)//":"//time_str(3:4)//":"//time_str(5:6)
        write(iu, '()')


        write(iu, '(A)') "Fit model:         " // trim(fitter_instance%model)
        write(iu, '(A,I0, " - ")', advance = 'no') "Likelihood:        ", fitter_instance%likelihood
        select case (fitter_instance%likelihood)
            case(1)  
                write(iu, '(A)') "Unbinned fit (chi2)"
            case(2)  
                write(iu, '(A)') "Binned fit (Poisson)"
            case(3)  
                write(iu, '(A)') "Binned fit (Gaussian)"
            case default
                write(iu, '(A)') "Unknown likelihood"
        end select

        write(iu, '(A,I0)') "Number of pars:    ", fitter_instance%npars

        write(iu, '(A)', advance = 'no') "Fit status:        "
        if(maxval(abs(fitter_instance%grad_p_min)) < fitter_instance%tol) then
            write(iu, '("CONVERGED after ", I0, " / " , I0 , " iterations.")') fitter_instance%niter, fitter_instance%max_iter
               
        else
            write(iu, '(A)', advance= 'no') "** NOT CONVERGED ** - gradient above tolerance for parameters: "
            do i = 1, fitter_instance%npars
                if(abs(fitter_instance%grad_p_min(i)) > fitter_instance%tol) then
                    write(iu, '(1x,A, 1x)', advance='no')  trim(fitter_instance%plabels(i))
                    
                end if
            end do
            
            write(iu, '(/,A, /)') "                 --> Please allow more iterations or fix parameter steps." 
        end if

        write(iu, '(A,F10.5,A,I0)') "FCN min value:   ", fitter_instance%min_FCN, "  / ndof ", fitter_instance%ndof
         
        write(iu, '(/,A,/)') "Parameter summary:"
        write(iu, '(A)') "----------------------------------------------------------------------------------"
        write(iu, '(A)') " Name     |  Fit Value   |  Uncertainty  |   Initial   |    Step    |    Grad "
        write(iu, '(A)') "----------------------------------------------------------------------------------"

        do j = 1, fitter_instance%npars
            write(iu, '(A9," |",F12.4,"  |",F12.4,"   |",F12.4," |",ES11.2," |",ES12.2)') trim(fitter_instance%plabels(j)), &
                fitter_instance%p_min(j), fitter_instance%p_err(j), fitter_instance%p0(j), &
                fitter_instance%step(j), fitter_instance%grad_p_min(j)
        end do

        write(iu, '(A)') "----------------------------------------------------------------------------------"



       
        write(iu, '(/,A)') "Covariance matrix  "
        write(iu, '(A)')   "------------------ "

        write(iu, '(A)', advance='no') "             "  
        do j = 1, fitter_instance%npars
            write(iu, '(A11,6x)', advance='no') adjustl(trim(fitter_instance%plabels(j)))
        end do
        write(iu, '(/)')

        do i = 1, fitter_instance%npars
            
            write(iu, '(A13)', advance='no') adjustl(trim(fitter_instance%plabels(i)))
            do j = 1, fitter_instance%npars
                if (j < i) then
                    write(iu, '(A17)', advance='no') " "  ! Blank lower triangle
                else
                    write(iu, '(ES15.3,2x)', advance='no') fitter_instance%cov_matrix(i,j)
                end if
            end do
            write(iu, '()')  ! End of row
        end do


        write(iu, '(/,A)') "Correlation matrix "
        write(iu, '(A)')   "------------------ "

        write(iu, '(A)', advance='no') "             "  
        do j = 1, fitter_instance%npars
            write(iu, '(A11,6x)', advance='no') adjustl(trim(fitter_instance%plabels(j)))
        end do
        write(iu, '(/)')

        do i = 1, fitter_instance%npars
            
            write(iu, '(A13)', advance='no') adjustl(trim(fitter_instance%plabels(i)))
            do j = 1, fitter_instance%npars
                if (j < i) then
                    write(iu, '(A17)', advance='no') " "  ! Blank lower triangle
                else
                    write(iu, '(F12.4, 5x)', advance='no') fitter_instance%cov_matrix(i,j) / &
                        (fitter_instance%p_err(i) * fitter_instance%p_err(j))
                end if
            end do
            write(iu, '()')  
        end do


        write(iu, '(/,A,/)') repeat('=', 90)


        if(present(filename)) then
            close(iu)
        end if
    

        
    end subroutine print_fit_results


    ! Export the histogram and function values to a file

    subroutine export_hist_and_fcn(fitter_instance, FCN, filename)
        
        interface
            function FCN(fitter_instance) result(val)
                import :: Fitter, dp
                type(Fitter) :: fitter_instance
                real(dp) :: val
            end function FCN
        end interface

        type(Fitter) :: fitter_instance
        character(len=*), intent(in) :: filename
        integer :: n, i, ipar, iu = 99, ncolumns = 0
        real(dp) ::  x, y, err_y

        if(fitter_instance%hist%nbins > 0) then
            n = fitter_instance%hist%nbins
        else
            n = size(fitter_instance%yvals)
        end if

        open(unit=iu, file=filename, status='replace')


        write(iu, '("# ", A13, 1x, A15, 1x)', advance = "no") "FCN value", "ndof"
        ncolumns = 2
        do ipar = 1, size(fitter_instance%p_min)
            write(iu, '(A15, 1x, A15, 1x)', advance = "no") &
                trim(fitter_instance%plabels(ipar)), "uncertainty"
            ncolumns = ncolumns + 2
        end do
        write(iu, '( 1x )')  
        write(iu, '("#", A) ') repeat('-', ncolumns * 16)

        fitter_instance%save_pdf = .true.          
        write(iu, '( F15.6, 1x, I15, 1x )', advance = "no") &
            FCN(fitter_instance), fitter_instance%ndof

        do ipar = 1, size(fitter_instance%p_min)
            write(iu, '( F15.6 , 1x, F15.6 , 1x)', advance = "no") &
                fitter_instance%p_min(ipar), fitter_instance%p_err(ipar)
        end do

        write(iu, '( /, 1x ,/ )')  

        

        write(iu, '("# ", A13, 1x, A15, 1x, A15, 1x, A15, 1x)') &
            "x", "y", "err_y", "f(x)"
        write(iu, '("#", A) ') repeat('-', ncolumns * 16)

        do i = 1, n

            if(fitter_instance%hist%nbins > 0) then
                x = fitter_instance%hist%bin_centers(i)
                y = fitter_instance%hist%bin_counts(i)
                err_y = fitter_instance%hist%bin_errors(i)
            else
                x = fitter_instance%xvals(i)
                y = fitter_instance%yvals(i)
                err_y = fitter_instance%err_yvals(i)
            end if
   
            write(iu, '(F15.6, 1x, F15.6, 1x, F15.6, 1x, F15.6, 1x )', advance = "no") &
                x, y, err_y, fitter_instance%pdf_vals(i) 
            
            do ipar = 5, ncolumns
                write(iu, '( F15.6, 1x )', advance = "no") &
                    0.0_dp  ! Placeholder for reading in python
            end do
            
            write(iu, '( 1x )')
            
        end do

        close(iu)

        write(*, '("Exported histogram and function values to ", A, " - to plot the results, use the following command:")')  trim(filename)
        write(*, '("  > python plot.py ", A )') trim(filename)

end subroutine export_hist_and_fcn


end module Minimizer
