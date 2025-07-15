!!!!!!!!!!!!!!!!!!!!!!!!! HistogramHandler.f08 !!!!!!!!!!!!!!!!!!!!!!!!!!
! This module handles histogram creation and manipulation

module HistogramHandler
    use utils

    type Histogram
        integer :: nbins
        
        real(dp), dimension(:), allocatable :: bin_edges
        real(dp), dimension(:), allocatable :: bin_centers
        real(dp), dimension(:), allocatable :: bin_widths

        integer :: nentries
        real(dp), dimension(:), allocatable :: bin_counts
        real(dp), dimension(:), allocatable :: bin_errors
        
        integer :: overflow_count 
        integer :: underflow_count
        character(len=100) :: label

    contains
        procedure :: fill_Histogram
        procedure :: print_Histogram
        
    end type Histogram
    

    
    interface Histogram
        module procedure init_Histogram
    end interface Histogram
    

    contains

        ! Initialize the histogram

        function init_Histogram(nbins, bin_min, bin_max, label) result(histogram_instance)
            integer, intent(in) :: nbins
            real(dp), intent(in) :: bin_min, bin_max
            character(len=*), intent(in), optional :: label
            integer :: i

            type(Histogram) :: histogram_instance

            if(bin_min >= bin_max) then
                write(*, '(/,A)') "Error: bin_min must be less than bin_max."
                stop
            end if
            if(nbins <= 0) then
                write(*, '(/,A)') "Error: nbins must be a positive integer."
                stop
            end if
            
            histogram_instance%nbins = nbins
            allocate(histogram_instance%bin_edges(nbins + 1))
            allocate(histogram_instance%bin_centers(nbins))
            allocate(histogram_instance%bin_widths(nbins))

            allocate(histogram_instance%bin_counts(nbins))
            allocate(histogram_instance%bin_errors(nbins))

            histogram_instance%nentries = 0
            histogram_instance%overflow_count = 0
            histogram_instance%underflow_count = 0

            if (present(label)) then
                histogram_instance%label = label
            else
                histogram_instance%label = "Histogram"
            end if
            histogram_instance%bin_edges = [(bin_min + i * (bin_max - bin_min) / real(nbins), i = 0, nbins)]
            histogram_instance%bin_centers = [(0.5_dp * (histogram_instance%bin_edges(i) + histogram_instance%bin_edges(i + 1)), i = 1, nbins)]
            histogram_instance%bin_widths = [(histogram_instance%bin_edges(i + 1) - histogram_instance%bin_edges(i), i = 1, nbins)]
            
            histogram_instance%bin_counts = 0.0_dp
            histogram_instance%bin_errors = 0.0_dp


        end function init_Histogram


        ! Fill the histogram with xvals

        subroutine fill_Histogram(this, n, xvals)
            class(Histogram) :: this
            integer, intent(in) :: n
            real(dp), dimension(n), intent(in) :: xvals
            integer :: i, ibin
            real(dp) :: bin_min, bin_max
            

            do i = 1, n
                if (xvals(i) < this%bin_edges(1)) then
                        this%underflow_count = this%underflow_count + 1
                else if (xvals(i) >= this%bin_edges(this%nbins + 1)) then
                        this%overflow_count = this%overflow_count + 1                
                else
                
                    do ibin = 1, this%nbins 
                        bin_min = this%bin_edges(ibin)
                        bin_max = this%bin_edges(ibin + 1)

                        if (xvals(i) >= bin_min .and. xvals(i) < bin_max) then
                            this%bin_counts(ibin) = this%bin_counts(ibin) + 1
                            this%nentries = this%nentries + 1
                            exit
                        end if
                    end do
                end if

            end do

            do ibin = 1, this%nbins
                ! Calculate errors for each bin (assuming Poisson statistics)
                if (this%bin_counts(ibin) > 0) then
                    this%bin_errors(ibin) = sqrt(this%bin_counts(ibin))
                else
                    this%bin_errors(ibin) = 0.0_dp
                end if
            end do

        end subroutine fill_Histogram

        ! Print the histogram contents

        subroutine print_Histogram(this)

            class(Histogram) :: this
            integer :: i, Nentries

        write(*, '(/,A)', advance = 'no') repeat('=', 35)
        write(*, '(A)', advance = 'no') " HISTOGRAM SUMMARY ="
        write(*, '(A,/)') repeat('=', 35)

       

        write(*, '(A, I0, A, F5.2, A, F5.2)') "  Number of bins:         ", this%nbins, "  from ", this%bin_edges(1), " to ", this%bin_edges(this%nbins + 1)
        write(*, '(A)') "  Histogram label:        " //  trim(this%label)

        write(*, '(A)', advance='no') "  Bin edges:              ["
        do i = 1, this%nbins + 1
            write(*, '(1x, F5.2)', advance='no') this%bin_edges(i)
            if (i < this%nbins + 1) write(*, '(A)', advance='no') ", "
            if (mod(i, 6) == 0 .and. i < this%nbins + 1) then
                print *, " "
                write(*, '(A)', advance='no') "                           "  
            end if
        end do
        print *, " ]"

        write(*, '(A, I0)') "  Total entries:          ", this%nentries
        write(*, '(A, I0)') "  Underflow count:        ", this%underflow_count
        write(*, '(A, I0)') "  Overflow count:         ", this%overflow_count
        print *, ""

        ! Print bin-by-bin content
        write(*, '(A)') "  Bin counts:"
        do i = 1, this%nbins
            write(*, '(A,I3,A,F6.2,A,F6.2,A,F5.0)') "    Bin ", i, ": [", this%bin_edges(i), ", ", this%bin_edges(i+1), "] -> ", this%bin_counts(i)
        end do

        write(*, '(/,A,/)') repeat('=', 90)



        end subroutine print_Histogram



end module HistogramHandler