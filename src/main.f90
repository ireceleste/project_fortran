program main

    use utils
    use InputOutput
    use HistogramHandler
    use GaussianGenerator
    use FitterModule
    use Minimizer
    use FunctionModule

    implicit none

    type(CardData) :: card    
    type(Fitter) :: fitter_instance

    ! ------------------Read input card and initialize the minimizer ------------------

    card = read_card()
    fitter_instance = init_fitter_from_card(card)

    ! ------------------ Perform minization and save results ------------------

    if(card%verbose_flag) then
         call gradient_minimize(fitter_instance, FCN, card%verbose_log)
    else
         call gradient_minimize(fitter_instance, FCN)
    end if
   
    call estimate_covariance(fitter_instance, FCN)

    call print_fit_results(fitter_instance)

    if(card%output_fit_summary_flag) then
        call print_fit_results(fitter_instance, card%output_fit_summary)
    end if

    if(card%output_fit_results_flag ) then
        call export_hist_and_fcn(fitter_instance, FCN, card%output_fit_results)
    end if

end program main
