module InputOutput

    use utils
    use FitterModule
    use HistogramHandler
    use GaussianGenerator
    use FunctionModule

    implicit none


    type CardData

        character(len=100) :: data_file, output_fit_summary, output_fit_results, verbose_log, model
        logical :: binned, input_data_flag, output_fit_summary_flag, output_fit_results_flag, verbose_flag, &
                    n_points_flag

        integer :: likelihood, Npars, max_iter, n_points
        real(dp), allocatable :: x_start(:), step(:)

    end type CardData


    contains

        function read_card() result(card)

            type(CardData) :: card

            
            character(len=100) :: cardname
            character(len=256) :: line
            character(len=32)  :: tokens(10)
            integer :: num_tokens, i, pos1, pos2, len_line
            
            integer :: iu = 10, ios, nargs
            logical :: exists
            logical :: likeflg, modelflg, start_params_flag, step_flag


            card%max_iter = 50000
            likeflg = .false.
            modelflg = .false.
            start_params_flag = .false.

            
            card%input_data_flag = .false. 
            card%output_fit_summary_flag = .false.
            card%output_fit_results_flag = .false.
            card%verbose_flag = .false.
            card%n_points_flag = .false.


            ! Check if the input card is provided

            nargs = command_argument_count()

            if (nargs == 0 .or. nargs >= 2) then
                write(*, '(/,A,/,A,/)') "Wrong input!", "Usage of the Fitter project:"
                write(*, '(A)') "  > make"
                write(*, '(A)') "  > bin/main path/to/inputcard"
                write(*, '(/,A,/)') "For help mode use:"
                write(*, '(A)') "  > make"
                write(*, '(A)') "  > bin/main -h  or  bin/main --help"
                write(*, '(/,A,/)') "For more details on the available models use:"
                write(*, '(A)') "  > make"
                write(*, '(A)') "  > bin/main -m  or  bin/main --models"
                write(*, '(/)')
                stop
            end if


            call get_command_argument(1,cardname)

            if (trim(cardname) == '--help' .or. trim(cardname) == '-h') then

                write(*, '(/,A,/)') repeat('=', 90)
                write(*, '(A)') "            GRADIENT DESCENT FITTER — HELP MODE"
                write(*, '(/,A,/)') repeat('=', 90)

                write(*, '(A)') "This program performs a numerical minimization of a fit function (chi2 or likelihood)"
                write(*, '(A)') "using the Gradient Descent method. It can operate on user-provided data or generate"
                write(*, '(A)') "synthetic datasets according to selected analytical models."
                write(*, '(/,A)') "The program supports:"
                write(*, '(A)') "  - Binned and unbinned fits"
                write(*, '(A)') "  - Gaussian or Poisson likelihood"
                write(*, '(A)') "  - Flexible model definitions (Gauss, Exp, Polynomial, etc.)"
                write(*, '(A)') "  - Numerical estimation of gradients and Hessians"
                write(*, '(A)') "  - Covariance matrix calculation using LAPACK (dgetrf, dgetri)"

                write(*, '(/,A)') "======================= INPUT FORMAT ======================="
                write(*, '(A)') "An input card must be provided as the first command-line argument:"
                write(*, '(A)') "  > bin/main path/to/inputcard"
                write(*, '(A)') "The input card is a two-column text file:"
                write(*, '(A)') "  [key]    [value]"
                write(*, '(A)') "Lines after the first empty line are ignored."

                write(*, '(/,A)') "Input options:"
                write(*, '(A)') "  likelihood              1 = unbinned (Gaussian), 2 = binned (Poisson), 3 = binned (Gaussian)"
                write(*, '(A)') "  model                   Gauss, Exp, Linear, Quadratic"
                write(*, '(A)') "  start_params            1.0 0.5 ! Initial parameters for the fit model (e.g., mean, stddev, slope, intercept, yield (optional))"
                write(*, '(A)') "  data_file               path/to/data.txt"
                write(*, '(A)') "  n_points                10000 ! Max number of points to generate (if no data_file is provided)"
                write(*, '(A)') "  output_fit_results      output/fit_results.txt"
                write(*, '(A)') "  fit_summary_file        output/fit_summary.txt"
                write(*, '(A)') "  verbose                 output/verbose:log.txt"

                write(*, '(/,A)') "======================= OUTPUT ======================="
                write(*, '(A)') "The program always prints to terminal:"
                write(*, '(A)') "  - Initial parameters, final fit values"
                write(*, '(A)') "  - Function value at minimum, number of degrees of freedom"
                write(*, '(A)') "  - Gradient, parameter errors, covariance and correlation matrix"
                write(*, '(A)') "Optional files to save the fit summary and results can be provided:"
                write(*, '(A)') "  - Fit summary (txt)"
                write(*, '(A)') "  - Fit results and input data for plotting (txt)"

                write(*, '(/,A,/)') repeat('=', 90)

                write(*, '(A)', advance='no') "Would you like to generate a template input card file (card.dat)? [y/n]: "
                read(*, '(A)') line

                if (trim(line) == 'y' .or. trim(line) == 'Y' .or. &
                    trim(line) == 'yes' .or. trim(line) == 'YES' .or. &
                    trim(line) == 'Yes') then


                    open(unit=iu, file='card.dat', status='replace', action='write')
                    write(iu,'(A)') "likelihood          2"
                    write(iu,'(A)') "model               Gauss"
                    write(iu,'(A)') "start_params        0.0 1.0 10000"
                    write(iu,'(A)') "data_file            "
                    write(iu,'(A)') "n_points             "
                    write(iu,'(A)') "max_iter            50000"
                    write(iu,'(A)') "step                1.0e-6 1.0e-6 10"
                    write(iu,'(A)') "output_fit_results  output/fit_results.txt"
                    write(iu,'(A)') "fit_summary_file    output/fit_summary.txt"
                    write(iu,'(A)') "verbose             output/verbose_log.txt"

                    write(iu,'(A)') ""
                    close(iu)


                    write(*, '(/,A,/)') "Template input card written to 'card.dat'."
                    write(*, '(A)') "Now you can modify it according to your needs and run the program with:"
                    write(*, '(A)') "  > bin/main card.dat"
                else
                    write(*, '(/,A,/)') "No file created. Proceed with your custom input card."
                end if

                stop

            else if (trim(cardname) == '--models' .or. trim(cardname) == '-m') then
                write(*, '(/,A,/)') repeat('=', 90)
                write(*, '(A)') "                     GRADIENT DESCENT FITTER — MODELS"
                write(*, '(/,A,/)') repeat('=', 90)

                write(*, '(A)') "Available models for UNBINNED fits (likelihood = 1):"
                write(*, '(A)') "  - Linear      (2 parameters: slope, intercept)"
                write(*, '(A)') "  - Quadratic   (3 parameters: a, b, c)"
                

                write(*, '(/,A)') "Available models for BINNED fits (likelihood = 2 = Poisson, or 3 = Gaussian):"
                write(*, '(A)') "  - Gaussian    (2 parameters: mean, stddev)"
                write(*, '(A)') "  - Exponential (1 parameter: lambda)"
                write(*, '(A)') "     > These models support an optional 3rd parameter: yield"
                write(*, '(A)') "       It can be used instead of the 'n_points' generation option."

                write(*, '(A)') "     > Parameter order for Gaussian: [mean, stddev, yield (optional)]"
                write(*, '(A)') "     > Parameter order for Exponential: [lambda, yield (optional)]"

                write(*, '(/,A,/)') repeat('=', 90)

                stop
            else
                write(*, '(/,A)') repeat('*', 90)
                write(*, '(A,/)') repeat('*', 90)
                write(*, '(A)') "            WELCOME TO THE GRADIENT DESCENT FITTER PROJECT"
                write(*, '(/,A)') repeat('*', 90)
                write(*, '(A,/)') repeat('*', 90)

            end if

            ! Read the input card file name from command line
            
            inquire(file=cardname, exist=exists)

            if (exists) then
                write(*, '(/,A,A,A,/)') "Reading '", trim(cardname), "' as input card..."
            else
                write(*, '(/,A,A,A)') "File '", trim(cardname), "' does not exist, aborting!"
                write(*, '(/)')
                stop
            end if


            open(newunit=iu, file=cardname, iostat = ios, iomsg=line, action='read')
            if (ios /= 0) then
                print *, line
                stop
            end if

        ! Loop through the input card and save until an empty line is found
        
        do
            read(iu, '(a)', iostat=ios) line

            if (ios == 0) then
                
                ! --- Read the line and save the key-value pairs ---
                
                tokens = ""
                num_tokens = 0
                len_line = len_trim(line)
                pos1 = 1

                do while (pos1 <= len_line .and. num_tokens < 10)

                    do while (pos1 <= len_line .and. line(pos1:pos1) == ' ')
                        pos1 = pos1 + 1
                    end do

                    if (pos1 > len_line) exit  ! Stop reading the file if the line is empty

                    ! Find next space 
                    pos2 = pos1
                    do while (pos2 <= len_line .and. line(pos2:pos2) /= ' ')
                        pos2 = pos2 + 1
                    end do

                    num_tokens = num_tokens + 1
                    tokens(num_tokens) = line(pos1:pos2-1)

                    pos1 = pos2
                end do


                if (num_tokens == 0) exit


                ! --- Process the key-value pairs from the input card ---

                if (tokens(1) == 'likelihood') then

                    if(num_tokens < 2) then
                        write(*, '(/,A)') "Input error: 'likelihood' option requires a value."
                        write(*, '(A)')   "The 'likelihood' option must be one of the following:"
                        write(*, '(A)')   "  1 --> Unbinned fit with Gaussian Likelihood"
                        write(*, '(A)')   "  2 --> Binned fit with Poisson Likelihood"
                        write(*, '(A)')   "  3 --> Binned fit with Gaussian Likelihood"
                        write(*, '(A,/)') "Aborting..."
                        stop
                    else if (num_tokens > 2) then
                        write(*, '(/,A)') "Input error: too many values provided for 'likelihood'."
                        write(*, '(A)')   "The 'likelihood' option must be one of the following:"
                        write(*, '(A)')   "  1 --> Unbinned fit with Gaussian Likelihood"
                        write(*, '(A)')   "  2 --> Binned fit with Poisson Likelihood"
                        write(*, '(A)')   "  3 --> Binned fit with Gaussian Likelihood"
                        write(*, '(A,/)') "Aborting..."
                        stop
                    end if

                    read(tokens(2), *, iostat=ios) card%likelihood

                    if (ios /= 0 .or. card%likelihood < 1 .or. card%likelihood > 3) then

                        write(*, '(/,A)') "Input error: invalid value provided for 'likelihood'."
                        write(*, '(A)')   "The 'likelihood' option must be one of the following:"
                        write(*, '(A)')   "  1 --> Unbinned fit with Gaussian Likelihood"
                        write(*, '(A)')   "  2 --> Binned fit with Poisson Likelihood"
                        write(*, '(A)')   "  3 --> Binned fit with Gaussian Likelihood"
                        write(*, '(A,/)') "Aborting..."
                        stop
                    else
                        likeflg = .true.
                        
                    end if
                end if ! End key == 'likelihood'

                if (tokens(1) == 'model') then

                    if (num_tokens < 2) then
                        write(*, '(/,A)') "Input error: 'model' option requires a value."
                        write(*, '(A)')   "For more details on available models use:"
                        write(*, '(A)')   "  > make"
                        write(*, '(A)')   "  > bin/main -m  or  bin/main --models"
                        write(*, '(A,/)') "Aborting..."
                        stop
                    else if (num_tokens > 2) then
                        write(*, '(/,A)') "Input error: too many values provided for 'model'."
                        write(*, '(A)')   "For more details on available models use:"
                        write(*, '(A)')   "  > make"
                        write(*, '(A)')   "  > bin/main -m  or  bin/main --models"
                        write(*, '(A,/)') "Aborting..."
                        stop
                    end if

                    if (trim(tokens(2)) == 'Linear' .or. trim(tokens(2)) == 'linear' .or. &
                        trim(tokens(2)) == 'poly1' .or. trim(tokens(2)) == 'LIN') then
                        card%model = 'linear'
                        modelflg = .true.
                    else if (trim(tokens(2)) == 'Quadratic' .or. trim(tokens(2)) == 'quadratic' .or. &
                                trim(tokens(2)) == 'poly2' .or. trim(tokens(2)) == 'QUAD') then
                        card%model = 'quadratic'
                        modelflg = .true.
                    else if (trim(tokens(2)) == 'Gauss' .or. trim(tokens(2)) == 'gauss' .or. &
                        trim(tokens(2)) == 'Gaussian' .or. trim(tokens(2)) == 'GAUSS') then
                        card%model = 'Gauss'
                        modelflg = .true.
                    else if (trim(tokens(2)) == 'Expo' .or. trim(tokens(2)) == 'expo' .or. trim(tokens(2)) == 'exp' .or. &
                             trim(tokens(2)) == 'Exponential' .or. trim(tokens(2)) == 'EXPO') then
                        card%model = 'expo'
                        modelflg = .true.
                    else
                        write(*, '(/,A)') "Input error: this model is not available."
                        write(*, '(A)')   "For more details on available models use:"
                        write(*, '(A)')   "  > make"
                        write(*, '(A)')   "  > bin/main -m  or  bin/main --models"
                        write(*, '(A,/)') "Aborting..."
                        stop
                    end if
                    
                end if ! End key == 'model'

                if(tokens(1) == 'start_params') then

                    if (num_tokens > 1) then
                        allocate(card%x_start(num_tokens - 1))
                        do i = 2, num_tokens
                            read(tokens(i), *, iostat=ios) card%x_start(i - 1)
                            if (ios /= 0) then
                                write(*, '(/,A, i0)') "Input error: invalid value provided for 'start_params' - position ", i-1
                                write(*, '(A)')   "Please provide valid numerical values."
                                write(*, '(A,/)') "Aborting..."
                                stop
                            end if
                            start_params_flag = .true.
                        end do
                    end if
                    
                end if ! End key == 'start_params'

                if(tokens(1) == 'step') then

                    if (num_tokens > 1) then
                        allocate(card%step(num_tokens - 1))
                        do i = 2, num_tokens
                            read(tokens(i), *, iostat=ios) card%step(i - 1)
                            if (ios /= 0) then
                                write(*, '(/,A, i0)') "Input error: invalid value provided for 'step' - position ", i-1
                                write(*, '(A)')   "Please provide valid numerical values."
                                write(*, '(A,/)') "Aborting..."
                                stop
                            end if
                            step_flag = .true.
                        end do
                    end if
                    
                end if ! End key == 'step'

                if (tokens(1) == 'data_file') then
                    select case (num_tokens)
                    case (1)
                        card%input_data_flag = .false.
                    case(2)
                        card%input_data_flag = .true.
                        card%data_file = trim(tokens(2))
                    case default
                        
                        write(*, '(/,A)') "Input error: invalid value provided for 'data_file'."
                        write(*, '(A)')   "Please provide a valid path to the data file."
                        write(*, '(A,/)') "Aborting..."
                        stop
                
                    end select
                    
                end if ! End key == 'data_file'

                if (tokens(1) == 'n_points') then
                    select case (num_tokens)
                    case (1)
                        card%n_points_flag = .false.
                    case(2)
                        card%n_points_flag = .true.
                        read(tokens(2), *, iostat=ios) card%n_points
                        
                        if (ios /= 0 .or. card%n_points <= 0) then
                            write(*, '(/,A)') "Input error: invalid value provided for 'n_points'."
                            write(*, '(A)')   "Please provide a positive integer value."
                            write(*, '(A,/)') "Aborting..."
                            stop
                        end if
                    case default
                        write(*, '(/,A)') "Input error: invalid value provided for 'n_points'."
                        write(*, '(A)')   "Please provide a valid positive integer value."
                        write(*, '(A,/)') "Aborting..."
                        stop
                
                    end select
                    
                end if ! End key == 'n_points'

                if (tokens(1) == 'fit_summary_file') then
                    select case (num_tokens)
                    case (1)
                        card%output_fit_summary_flag = .false.
                    case(2)
                        card%output_fit_summary_flag = .true.
                        card%output_fit_summary = trim(tokens(2))
                    case default
                        
                        write(*, '(/,A)') "Input error: invalid value provided for 'fit_summary_file'."
                        write(*, '(A)')   "Please provide a valid path to the data file."
                        write(*, '(A,/)') "Aborting..."
                        stop
                
                    end select
                    
                end if ! End key == 'fit_summary_file'

                if (tokens(1) == 'output_fit_results') then
                    select case (num_tokens)
                    case (1)
                        card%output_fit_results_flag = .false.
                    case(2)
                        card%output_fit_results_flag = .true.
                        card%output_fit_results = trim(tokens(2))
                    case default
                        
                        write(*, '(/,A)') "Input error: invalid value provided for 'output_fit_results'."
                        write(*, '(A)')   "Please provide a valid path to the data file."
                        write(*, '(A,/)') "Aborting..."
                        stop
                
                    end select
                    
                end if ! End key == 'output_fit_results'
                
                if (tokens(1) == 'verbose') then
                    select case (num_tokens)
                    case (1)
                        card%verbose_flag = .false.
                    case(2)
                        card%verbose_flag = .true.
                        card%verbose_log = trim(tokens(2))
                    case default
                        
                        write(*, '(/,A)') "Input error: invalid value provided for 'verbose'."
                        write(*, '(A)')   "Please provide a valid path to the data file."
                        write(*, '(A,/)') "Aborting..."
                        stop
                
                    end select
                    
                end if ! End key == 'verbose'

                if(tokens(1) == 'max_iter' .and. num_tokens == 2) then
                    
                    read(tokens(2), *, iostat=ios) card%max_iter
                    if (ios /= 0 .or. card%max_iter <= 0) then
                        write(*, '(/,A)') "Input error: invalid value provided for 'max_iter'."
                        write(*, '(A)')   "Please provide a positive integer value."
                        write(*, '(A,/)') "Aborting..."
                        stop
                    end if
                                        
                end if ! End key == 'max_iter'


            else
                exit
            end if
        end do  ! End of input card reading loop



        ! Infer model properties

        if (likeflg .and. modelflg) then
            
            card%binned = (card%likelihood == 2 .or. card%likelihood == 3)

            if(card%binned) then
                if (card%model  == "Gauss") then
                    card%Npars = 3
                   
                else if (card%model  == "expo") then
                    card%Npars = 2
                    
                else
                    write(*, '(/,A)') "Input error: BINNED model not supported."
                    stop 
                end if
            else
                if (card%model  == "linear") then
                    card%Npars = 2
                    
                else if (card%model  == "quadratic") then
                    card%Npars = 3
                else
                    write(*, '(/,A)') "Input error: UNBINNED model not supported."
                    stop 
                end if
            end if
        else
            write(*, '(/,A)') "Input error: both likelihood and model must be specified."
            stop 
        end if

        ! Check if start parameters are provided in the correct format
        
        if(start_params_flag) then

            ! Stop if the number of start parameters does not match the expected number
            if(.not.card%binned .and. size(card%x_start) /= card%Npars) then
                write(*, '(/,A)') "Input error: number of start parameters does not match the expected number for unbinned fits."
                write(*, '("Expected ", I0, " parameters but got ", I0)')  card%Npars, size(card%x_start)
                stop
            else if(card%binned .and. (size(card%x_start) /= card%Npars - 1 .and. size(card%x_start) /= card%Npars )) then
                write(*, '(/,A)') "Input error: number of start parameters does not match the expected number for binned fit."
                write(*, '("Expected ", I0, " parameters but got ", I0)')  card%Npars - 1, size(card%x_start)
                stop
            
            end if
            
        else
            write(*, '(/,A)') "No start parameters provided, using default values."
            
            if(.not.card%n_points_flag) then
                if(card%binned) then
                    card%n_points = 10000  ! Default number of points for binned data generation
                else
                    card%n_points = 50  ! Default number of points for unbinned data generation
                end if
            end if

            allocate(card%x_start(card%Npars))
            
            select case (card%model)
            case ("Gauss")  ! Gauss

                card%x_start(1) = 1.0_dp  ! mean
                card%x_start(2) = 0.5_dp  ! stddev
                
            case ("expo")  ! Exponential
                card%x_start(1) = 1.0_dp  ! lambda
            case ("linear")  ! Linear
                card%x_start = [1.0_dp, 0.5_dp]  ! slope, intercept
            case ("quadratic")  ! Quadratic
                card%x_start = [1.0_dp, 0.5_dp, 0.1_dp]  ! a, b, c
            case default
                write(*, '(/,A)') "Input error: invalid model."
                stop
            end select
            
            
        end if

        ! Set correct number of points for binned fits if not given in the input card

        if(.not. card%n_points_flag ) then

            if(card%binned) then

                if( .not. start_params_flag  ) then
                    card%n_points = 10000  ! Default number of points for binned data generation
                    card%x_start(card%Npars) = real(card%n_points, dp)  
                
                else if(size(card%x_start) == card%Npars - 1) then
                    card %n_points = 10000  ! Default number of points for binned data generation
                    call append_element(card%x_start, real(card%n_points, dp) )

                else if(size(card%x_start) == card%Npars) then
                    card%n_points = int(card%x_start(card%Npars))  ! Set n_points from the yield parameter
                end if

            else
                card%n_points = 50  ! Default number of points for unbinned data generation
            end if
            
        end if

        ! Set the step if not provided in the input card

        if (.not. step_flag) then
            allocate(card%step(card%Npars))
            card%step = 1.0e-6_dp  ! Default step size for all parameters

            select case (card%model)
            case ("Gauss")  ! Gauss
                card%step(3) = 10.0_dp  ! yield step to ensure better convergence
            case ("expo")  ! Exponential
                card%step(2) = 10.0_dp  ! yield step to ensure better convergence
            case ("linear")  ! Linear
                ! No specific step size needed for linear model
            case ("quadratic")  ! Quadratic
                card%step(1) = 1.0e-9_dp  ! a step to ensure better convergence
            end select
            
        end if

    
        close(unit=iu)       


        ! Print a recap of the input card

        write(*, '(/,A)', advance = 'no') repeat('=', 34)
        write(*, '(A)', advance = 'no') " RECAP OF INPUT CARD "
        write(*, '(A,/)') repeat('=', 34)


        write(*, '(A,1X,A)', advance   ='no') "  Likelihood type:       "
        select case (card%likelihood)
        case (1)
            write(*, '(A)') "Unbinned (Gaussian)"
        case (2)
            write(*, '(A)') "Binned (Poisson)"
        case (3)
            write(*, '(A)') "Binned (Gaussian)"
        end select



        write(*, '(A,1X,A)') "  Model:                 ", trim(card%model)

        if(start_params_flag) then
            write(*, '("  Start parameters:       [")', advance='no')
            do i = 1, size(card%x_start)
                write(*, '(F8.3)', advance='no') card%x_start(i)
                if (i < size(card%x_start)) write(*, '(" ")', advance='no')
            end do
            write(*, '("]")')
        else
            write(*, '(A)') "  Start parameters:       [not set] -> will use default values:"
            write(*, '("                          [")', advance='no')
            do i = 1, size(card%x_start)
                write(*, '(F8.3)', advance='no') card%x_start(i)
                if (i < size(card%x_start)) write(*, '(" ")', advance='no')
            end do
            write(*, '("]")')
        end if



        write(*, '(A,1X,I0)')            "  Max iterations:        ", card%max_iter

        
        if (card%input_data_flag) then
            write(*, '(A,1X,A)') "  Data file:             ", trim(card%data_file)
        else
            write(*, '(A, I0)') "  Data file:              [none] -> will generate 'n_points' = ",  card%n_points 
        end if


        if (card%output_fit_summary_flag) then
            write(*, '(A,1X,A)') "  Fit summary file:      ", trim(card%output_fit_summary)
        else
            write(*, '(A)') "  Fit summary file:       [not set] -> will not save fit summary"
        end if

        if (card%output_fit_results_flag) then
            write(*, '(A,1X,A)') "  Fit results file:      ", trim(card%output_fit_results)
        else
            write(*, '(A)') "  Fit results file:       [not set] -> will not save fit results"
        end if

        if (card%verbose_flag) then
            write(*, '(A,1X,A)') "  Verbose log file:      ", trim(card%verbose_log)
        else
            write(*, '(A)') "  Verbose log file:       [not set] -> will not save verbose log"
        end if

        write(*, '(/,A,/)') repeat('=', 90)



    end function read_card

    
    function init_fitter_from_card(card) result(fitter_instance)

        type(CardData) :: card
        type(Fitter) :: fitter_instance
        
        integer ::  n, i, ios, iu = 11, nbins_hist
        real(dp) :: rng, xmin, xmax
        character(len=256) :: line

        real(dp), dimension(:), allocatable :: xvals, yvals, err_yvals
        real(dp), dimension(:), allocatable :: lrate
        character(len=100), dimension(:), allocatable :: param_names
        character(len=100) :: hist_name
        type(Histogram) :: hist
        logical :: exists

        ! Set step, learning rate, and parameter names based on the FCN index


        allocate(lrate(card%Npars), param_names(card%Npars))
        
        select case (card%model)
        case ("Gauss")  
            param_names(1) = "mean"
            param_names(2) = "stddev"
            param_names(3) = "yield"
        case ("expo")
            param_names(1) = "lambda"
            param_names(2) = "yield"
        case ("linear") 
            param_names(1) = "slope"
            param_names(2) = "intercept"
        case ("quadratic")  
            param_names(1) = "a"
            param_names(2) = "b"
            param_names(3) = "c"
        case default
            write(*, '(/,A)') "Input error: invalid card model."
            stop
        end select

        lrate = card%step  ! Use the step sizes as learning rates


        ! Read or generate data 

        if(card%binned) then

            if(card%input_data_flag) then
            

                inquire(file=card%data_file, exist=exists)
                if (.not. exists) then
                    write(*, '(/,A,A)') "Input error: data file '", trim(card%data_file), "' does not exist."
                    stop
                end if

                open(newunit=iu, file=card%data_file, iostat = ios, iomsg=line, action='read', status='old')

                if (ios /= 0) then
                    write(*, '(/,A,A)') "Input error: cannot open data file '", trim(card%data_file), "'."
                    stop
                end if
                
                n = 0
                do
                    read(iu, *, iostat=ios) xvals(i)
                    if (ios /= 0) exit
                    n = n + 1
                end do
                write(*, '(A,I0)') "Read ", n, " data points from file '", trim(card%data_file), "'."
                hist_name = "Histogram of " // trim(card%data_file) // " data"
                
            else
                n = card%n_points
                allocate(xvals(n))
                call generate_dataset_pdf(n, xvals, pdf = card%model, pars = card%x_start(1:card%Npars-1)) 

                hist_name = "Histogram of generated " // trim(card%model) // " data"

            end if

            card%x_start(card%Npars) = real(n, dp)  ! Update yield as the number of points to ensure convergence

            if(card%model == "Gauss") then
                xmin = card%x_start(1) - 3.0_dp * card%x_start(2)
                xmax = card%x_start(1) + 3.0_dp * card%x_start(2)
            else if (card%model == "expo") then
                xmin = 0.0_dp
                xmax = 3.0_dp * card%x_start(1)
            end if

            nbins_hist = 50

            if(card%n_points < 1000) then
                nbins_hist = 10  
            end if

             
            hist = Histogram(nbins_hist, xmin, xmax,  hist_name)
            call hist%fill_Histogram(n, xvals)
            call hist%print_Histogram()

            fitter_instance = init_Fitter(card%likelihood, card%model, card%Npars, &
                card%x_start, param_names, &
                hist=hist , &
                max_iter=card%max_iter, step=card%step, lrate=lrate)

        
        else ! Unbinned fit
            if(card%input_data_flag) then
                inquire(file=card%data_file, exist=exists)
                if (.not. exists) then
                    write(*, '(/,A,A)') "Input error: data file '", trim(card%data_file), "' does not exist."
                    stop
                end if

                open(newunit=iu, file=card%data_file, iostat = ios, iomsg=line, action='read', status='old')

                if (ios /= 0) then
                    write(*, '(/,A,A)') "Input error: cannot open data file '", trim(card%data_file), "'."
                    stop
                end if
                
                n = 0
                do
                    read(iu, *, iostat=ios) xvals(n), yvals(n), err_yvals(n)
                    if (ios /= 0) exit
                    n = n + 1
                end do
                write(*, '(A,I0)') "Read ", n, " data points from file '", trim(card%data_file), "'."
                
            else
                n = card%n_points
                allocate(xvals(n), yvals(n), err_yvals(n))
                do i = 1, n
                    call random_number(rng)
                    xvals(i) = -10.0_dp + rng * 20.0_dp
                    call random_number(rng)
                    err_yvals(i) = 0.1_dp + rng * 0.5_dp  ! Random error between 0.1 and 0.6
                end do
        
                call generate_dataset_model(n, yvals, card%x_start, model = card%model, &
                     xvals = xvals, yerrs = err_yvals)
            end if

            fitter_instance = init_Fitter(card%likelihood, card%model, card%Npars, &
                card%x_start, param_names, &
                xvals = xvals, yvals=yvals, err_yvals=err_yvals, &
                max_iter=card%max_iter, step=card%step, lrate=lrate)
        end if

       

    end function init_fitter_from_card


    
end module InputOutput