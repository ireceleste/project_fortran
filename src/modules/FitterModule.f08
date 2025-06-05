

module FitterModule
    use utils, only: dp  
    use HistogramHandler
    implicit none


    type Fitter

        ! -------- Fitter object to hold parameters and yvals for fitting --------

        ! Parameters
        ! index: identifier for the fit strategy: 
        !         1: Gaussian fit to histogram
        ! npars: number of parameters to fit
        ! p0: initial parameter values
        ! optional, plabels: labels for the parameters
        ! optional, hist: Histogram object to hold histogram yvals
        ! optional, yvals: yvals points to fit
        ! optional, err_yvals: errors associated with the yvals points
        ! optional, xvals: x values corresponding to the yvals points
        ! optional, save_pdf: flag to save PDF values
        ! tol: tolerance for convergence
        ! max_iter: maximum number of iterations for fitting
        ! step: step size for parameter updates
        ! lrate: learning rate for optimization 


        integer :: likelihood
        character(len=100) :: model
        
        integer :: npars 
        real(dp), dimension(:), allocatable :: p0
        character(len=100) , dimension(:), allocatable :: plabels

        type(Histogram) :: hist
        real(dp), dimension(:), allocatable :: yvals, err_yvals, xvals
        real(dp), dimension(:), allocatable :: pdf_vals
        logical :: save_pdf 
        
        
        real(dp) :: tol
        integer :: max_iter
        real(dp), dimension(:), allocatable :: step
        real(dp), dimension(:), allocatable :: lrate

        
        ! Results
        ! p_min: minimum parameters found during fitting
        ! grad_p_min: gradient of the fitting function at the minimum parameters
        ! cov_matrix: covariance matrix of the fitted parameters
        ! p_err: errors associated with the fitted parameters
        ! min_FCN: minimum value of the fitting function
        ! ndof: number of degrees of freedom in the fit
        ! niter: number of iterations taken to converge

        real(dp), dimension(:), allocatable :: p_min
        real(dp), dimension(:), allocatable :: grad_p_min
        real(dp), dimension(:,:), allocatable :: cov_matrix
        real(dp), dimension(:), allocatable :: p_err

        real(dp) :: min_FCN 
        integer :: ndof
        integer :: niter

    end type Fitter


    contains

    ! Initialize the Fitter object 

    function init_Fitter( likelihood, model, npars, p0, plabels, hist, yvals, err_yvals, xvals, tol, max_iter, step, lrate) result(fitter_instance)
       
        implicit none

        type(Fitter) :: fitter_instance
        integer :: i

        ! Necessary input parameters

        integer, intent(in) :: likelihood, npars
        character(len=*), intent(in) :: model
        real(dp), dimension(:), intent(in) :: p0


        ! Optional parameters
        character(len=*), dimension(:), intent(in), optional :: plabels
        real(dp), dimension(:), intent(in), optional :: step
        real(dp), dimension(:), intent(in), optional :: lrate
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: max_iter

        type(Histogram), intent(in), optional :: hist
        real(dp), dimension(:), intent(in), optional :: yvals, err_yvals, xvals



        ! --- Input parameters initialisation ---

        fitter_instance%likelihood = likelihood
        fitter_instance%model = model
        fitter_instance%npars = npars

        ! Allocate and initialize parameters
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

        ! yvals initialization
        if (present(yvals)) then
            allocate(fitter_instance%yvals(size(yvals)))
            fitter_instance%yvals = yvals
        else
            allocate(fitter_instance%yvals(0))  
        end if
        
        ! Error yvals initialization
        if (present(err_yvals)) then
            allocate(fitter_instance%err_yvals(size(err_yvals)))
            fitter_instance%err_yvals = err_yvals
        else
            allocate(fitter_instance%err_yvals(size(yvals))) 
        end if

        ! X values initialization
        if (present(xvals)) then
            allocate(fitter_instance%xvals(size(xvals)))
            fitter_instance%xvals = xvals
        else
            allocate(fitter_instance%xvals(size(yvals))) 
        end if


        ! Model values initialization
        if(present(xvals)) then

            allocate(fitter_instance%pdf_vals(size(xvals)))
            fitter_instance%pdf_vals = 0.0_dp  

        else if(present(hist) ) then
            allocate(fitter_instance%pdf_vals(hist%nbins))
            fitter_instance%pdf_vals = 0.0_dp  
        else
            allocate(fitter_instance%pdf_vals(size(yvals)))  
        end if

        ! Save PDF flag initialization
        fitter_instance%save_pdf = .false. 


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
            fitter_instance%step = 1.0e-6_dp  ! Default step size
        end if

        ! Learning rate initialization
        if (present(lrate)) then
            allocate(fitter_instance%lrate(size(lrate)))
            fitter_instance%lrate = lrate
        else
            allocate(fitter_instance%lrate(fitter_instance%npars))
            fitter_instance%lrate = 1.0e-6_dp   ! Default learning rate
        end if

        ! Ensure the learning rate and the steps are positive

        do i = 1, size(fitter_instance%lrate)
            if (fitter_instance%lrate(i) <= 0.0_dp) then
                write(*, '( "Warning: Learning rate must be positive: Index " , I0, " is " , F10.5, " ---> Setting to default value of 0.01.")') i , fitter_instance%lrate(i)
                fitter_instance%lrate(i) = 1.0e-6_dp 
            end if
            if (fitter_instance%step(i) <= 0.0_dp) then
                write(*, '( "Warning: Step size must be positive: Index " , I0, " is " , F10.5, " ---> Setting to default value of 0.01.")') i , fitter_instance%step(i)
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
        
        ! Initialize the covariance matrix
        allocate(fitter_instance%cov_matrix(fitter_instance%npars, fitter_instance%npars))
        fitter_instance%cov_matrix = 0.0_dp

        ! Initialize the parameter errors
        allocate(fitter_instance%p_err(fitter_instance%npars))
        fitter_instance%p_err = 0.0_dp


        ! Initialize the minimum function value
        fitter_instance%min_FCN = 0.0_dp

        ! Initialize degrees of freedom
        if(present(yvals)) then
            fitter_instance%ndof = size(yvals) - npars
        else if(present(hist)) then
            fitter_instance%ndof = hist%nbins - npars
        else
            fitter_instance%ndof = 0  ! Default case if no yvals or histogram is provided
        end if

        ! Initialize the number of iterations
        fitter_instance%niter = 0


    end function init_Fitter

end module FitterModule