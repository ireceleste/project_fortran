!!!!!!!!!!!!!!!!!!!!!!!!! FitterModule.f08 !!!!!!!!!!!!!!!!!!!!!!!!!!
! This module defines the Fitter object and contains functions to initialize it.
! A Fitter contains all the parameters and results of a fit.

module FitterModule
    use utils, only: dp  
    use HistogramHandler
    implicit none


    type Fitter

        ! -------- Fitter object that contains fit settings and results --------

        ! -- Input parameters, for fit settings --
            ! likelihood: identifier for the likelihood function to use
            !            1: Unbinned fit
            !            2: Binned fit with Poisson likelihood
            !            3: Binned fit with Gaussian likelihood
            ! model: string representing the model to fit
            ! npars: number of parameters to fit
            ! p0: initial parameter values

            ! optional, plabels: labels for the parameters
            ! optional, hist: Histogram object to hold histogram yvals
            ! optional, xvals: x values corresponding to the yvals points
            ! optional, yvals: yvals points to fit
            ! optional, err_yvals: errors associated with the yvals points
            ! optional, tol: tolerance for convergence
            ! optional, max_iter: maximum number of iterations for fitting
            ! optional, step: step size for parameter updates
            ! optional, lrate: learning rate for optimization 


        integer :: likelihood
        character(len=100) :: model
        integer :: npars 
        real(dp), dimension(:), allocatable :: p0

        character(len=100) , dimension(:), allocatable :: plabels
        type(Histogram) :: hist
        real(dp), dimension(:), allocatable :: xvals, yvals, err_yvals
        real(dp) :: tol
        integer :: max_iter
        real(dp), dimension(:), allocatable :: step, lrate

        
        ! -- Output parameters, for fit results --
            ! p_min: minimum parameters found during fitting
            ! grad_p_min: gradient of the fitting function at the minimum parameters
            ! p_err: errors associated with the fitted parameters
            ! cov_matrix: covariance matrix of the fitted parameters
            ! min_FCN: minimum value of the fitting function
            ! ndof: number of degrees of freedom in the fit
            ! niter: number of iterations taken to converge
            ! save_pdf: flag to save PDF values (set to true when the fit is performed)
            ! pdf_vals: PDF values corresponding to the yvals points

        real(dp), dimension(:), allocatable :: p_min, grad_p_min, p_err
        real(dp), dimension(:,:), allocatable :: cov_matrix
        real(dp) :: min_FCN 
        integer :: ndof, niter
        logical :: save_pdf 
        real(dp), dimension(:), allocatable :: pdf_vals


    end type Fitter


    contains

    ! Initialize the Fitter object 

    function init_Fitter( likelihood, model, npars, p0, plabels, hist, yvals, err_yvals, xvals, tol, max_iter, step, lrate) result(fitter_instance)
       
        implicit none

        type(Fitter) :: fitter_instance
        integer :: i

        ! Necessary input parameters
        integer, intent(in) :: likelihood
        character(len=*), intent(in) :: model
        integer, intent(in) :: npars
        real(dp), dimension(:), intent(in) :: p0


        ! Optional input parameters
        character(len=*), dimension(:), intent(in), optional :: plabels
         type(Histogram), intent(in), optional :: hist
        real(dp), dimension(:), intent(in), optional :: xvals, yvals, err_yvals
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: max_iter
        real(dp), dimension(:), intent(in), optional :: step, lrate

       
        ! --- Input parameters initialisation ---

        fitter_instance%likelihood = likelihood
        fitter_instance%model = model
        fitter_instance%npars = npars
        allocate(fitter_instance%p0(fitter_instance%npars))
        fitter_instance%p0 = p0

        ! Allocate and initialize parameter labels
        if(present(plabels)) then
            allocate(fitter_instance%plabels(size(plabels)))
            fitter_instance%plabels = plabels
        else
            allocate(fitter_instance%plabels(fitter_instance%npars))
            do i = 1, npars
                write(fitter_instance%plabels(i), '(A,I0)') 'Parameter ', i
            end do        
        end if


        ! Histogram initialization
        if (present(hist)) then
            fitter_instance%hist = hist
        else
            fitter_instance%hist = Histogram(0, 1.0_dp, 1.0_dp, "Null Histogram")
        end if

        ! xvals and yvals initialization
        if (present(yvals) .and. present(xvals)) then
            if(size(xvals) == size(yvals)) then
                allocate(fitter_instance%yvals(size(yvals)))
                fitter_instance%yvals = yvals
                allocate(fitter_instance%xvals(size(xvals)))
                fitter_instance%xvals = xvals
            else
                write(*,*) "Error: yvals and xvals must have the same size."
                stop
            end if
           
        else
            allocate(fitter_instance%yvals(0))  
            allocate(fitter_instance%xvals(0))
        end if
        
        ! Error yvals initialization
        if (present(err_yvals)) then
            allocate(fitter_instance%err_yvals(size(err_yvals)))
            fitter_instance%err_yvals = err_yvals
        else if (present(yvals) .and. present(xvals)) then
            allocate(fitter_instance%err_yvals(size(yvals))) 
            fitter_instance%err_yvals = 1.0_dp
        else
            allocate(fitter_instance%err_yvals(0))  
        end if

        

        ! Tolerance initialization
        if (present(tol)) then
            fitter_instance%tol = tol
        else
            fitter_instance%tol = 1.0e-6_dp  
        end if

        ! Maximum iterations initialization
        if (present(max_iter)) then
            fitter_instance%max_iter = max_iter
        else
            fitter_instance%max_iter = 1000  
        end if

        ! Step initialization
        if (present(step)) then
            allocate(fitter_instance%step(size(step)))
            fitter_instance%step = step
        else
            allocate(fitter_instance%step(fitter_instance%npars))
            fitter_instance%step = 1.0e-6_dp  
        end if

        ! Learning rate initialization
        if (present(lrate)) then
            allocate(fitter_instance%lrate(size(lrate)))
            fitter_instance%lrate = lrate
        else
            allocate(fitter_instance%lrate(fitter_instance%npars))
            fitter_instance%lrate = 1.0e-6_dp   
        end if

        ! Ensure the learning rate and the steps are positive

        do i = 1, size(fitter_instance%lrate)
            if (fitter_instance%lrate(i) <= 0.0_dp) then
                write(*, '( "Warning: Learning rate must be positive: Index " , I0, " is " , F10.5, " ---> Setting to default value of 1e-6.")') i , fitter_instance%lrate(i)
                fitter_instance%lrate(i) = 1.0e-6_dp 
            end if
            if (fitter_instance%step(i) <= 0.0_dp) then
                write(*, '( "Warning: Step size must be positive: Index " , I0, " is " , F10.5, " ---> Setting to default value of 1e-6.")') i , fitter_instance%step(i)
                fitter_instance%step(i) = 1.0e-6_dp 
            end if
        end do



        ! ---- Output parameters initialization ----

        ! Initialize the minimum parameters
        allocate(fitter_instance%p_min(size(p0)))
        fitter_instance%p_min = p0

        ! Initialize the gradient at the minimum parameters
        allocate(fitter_instance%grad_p_min(fitter_instance%npars))
        fitter_instance%grad_p_min = 0.0_dp
        
        ! Initialize the parameter errors
        allocate(fitter_instance%p_err(fitter_instance%npars))
        fitter_instance%p_err = 0.0_dp

        ! Initialize the covariance matrix
        allocate(fitter_instance%cov_matrix(fitter_instance%npars, fitter_instance%npars))
        fitter_instance%cov_matrix = 0.0_dp

        ! Initialize the minimum function value
        fitter_instance%min_FCN = 0.0_dp

        ! Initialize degrees of freedom
        if(present(yvals) .and. present(xvals)) then
            fitter_instance%ndof = size(yvals) - npars
        else if(present(hist)) then
            fitter_instance%ndof = hist%nbins - npars
        else
            fitter_instance%ndof = 0  ! Default case if no yvals or histogram is provided
        end if

        if (fitter_instance%ndof < 0) then
            write(*,'("Warning: Negative number of degrees of freedom!")')
            stop
        end if

        ! Initialize the number of iterations
        fitter_instance%niter = 0
        
        ! Save PDF flag initialization. Will be set to true when the fit is performed   
        fitter_instance%save_pdf = .false. 

        ! Saved pdf values initialization
        if(present(xvals) .and. present(yvals)) then
            allocate(fitter_instance%pdf_vals(size(xvals)))
            fitter_instance%pdf_vals = 0.0_dp  
        else if(present(hist) ) then
            allocate(fitter_instance%pdf_vals(hist%nbins))
            fitter_instance%pdf_vals = 0.0_dp  
        else
            allocate(fitter_instance%pdf_vals(0))
        end if


    end function init_Fitter

end module FitterModule