module utils
    integer, parameter :: dp = kind(1.0d0)
    real(dp) :: pi = 4.0_dp * atan(1.0_dp)
    real(dp) :: e = exp(1.0_dp)
    
    interface xgetrf
        module procedure dgetrf_interface
    end interface

    interface xgetri
        module procedure dgetri_interface
    end interface

    contains    

        ! Subroutine to append an element to an allocatable array
        
        subroutine append_element(array, value)
            implicit none
            real(dp), allocatable, intent(inout) :: array(:)
            real(dp), intent(in) :: value
            real(dp), allocatable :: temp(:)
            integer :: n

            n = size(array)
            allocate(temp(n))
            temp = array

            deallocate(array)
            allocate(array(n+1))
            array(1:n) = temp
            array(n+1) = value
        end subroutine append_element


        ! Function to check if two real numbers are close to each other

        logical function closeto(x, y, relative, absolute )

            real(dp) , intent(in):: x,y
            real , intent(in) , optional :: relative, absolute
            real :: th_relative, th_absolute, threshold

            ! Assign the thresholds using the default values or the input values
            if(present(relative)) then 
                th_relative = relative
            else
                th_relative = 1e-10_dp 
            end if

            if(present(absolute)) then 
                th_absolute = absolute
            else
                th_absolute = 1e-10_dp
            end if

            ! Now choose the global threshold and test both absolute and relative differences

            threshold = min(th_absolute, th_relative * max(abs(x), abs(y)))
            if(max(abs(x), abs(y)) < th_absolute) then
                threshold = th_absolute
            end if

            closeto = abs(x - y) .le. threshold

        end function closeto

        ! Interface wrapper for DGETRF 
        subroutine dgetrf_interface(m, n, a, lda, ipiv, info)
            integer, intent(in) :: m, n, lda
            integer, intent(out) :: ipiv(*)
            integer, intent(out) :: info
            real(dp), intent(inout) :: a(lda, *)
            external :: dgetrf
            call dgetrf(m, n, a, lda, ipiv, info)
        end subroutine dgetrf_interface

        ! Interface wrapper for DGETRI 
        subroutine dgetri_interface(n, a, lda, ipiv, work, lwork, info)
            integer, intent(in) :: n, lda, lwork
            integer, intent(in) :: ipiv(*)
            integer, intent(out) :: info
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(inout) :: work(*)
            external :: dgetri
            call dgetri(n, a, lda, ipiv, work, lwork, info)
        end subroutine dgetri_interface


        ! Procedure to invert a matrix using LAPACK's dgetrf and dgetri routines
        
        subroutine invert_matrix(A, Ainv)
        

            real(dp), intent(in)  :: A(:,:)
            real(dp), intent(out) :: Ainv(:,:)
            integer :: n, info, lwork
            integer, allocatable :: ipiv(:)
            real(dp), allocatable :: work(:), identity(:,:)
            logical :: success
            real(dp) :: exp


            n = size(A, 1)
            Ainv = A
            allocate(ipiv(n))

            ! First, perform LU decomposition
            call xgetrf(n, n, Ainv, n, ipiv, info) 

            if (info /= 0) then
                print*, "################### Error in LU decomposition ###################"
                print *, " LU decomposition failed. INFO =", info
                print*, "#################################################################"
                stop
            end if


            lwork = -1
            allocate(work(1))
            call xgetri(n, Ainv, n, ipiv, work, lwork, info)
            lwork = int(work(1))  
            deallocate(work)
            allocate(work(lwork))

            ! Perform matrix inversion
            call dgetri(n, Ainv, n, ipiv, work, lwork, info)

            if (info /= 0) then
                print *, "################### Error in matrix inversion ###################"
                print *, " Matrix inversion failed. INFO =", info
                print *, "#################################################################"
                stop
            end if

            deallocate(ipiv, work)

            ! check if Ainv is the inverse of A

            allocate(identity(n,n))
            identity = matmul(A, Ainv)


            do i = 1, n
                do j = 1, n
                   
                    if (i == j) then
                        exp = 1.0_dp  
                    else
                        exp = 0.0_dp
                    end if

                    if ( closeto(identity(i,j), exp ) ) then
                        success = .true.
                    else
                        success = .false.
                        print *, "################### Inversion check failed ###################"
                        write(*, '(A,I0,A,I0,A,F20.10)') " Identity matrix element (", i, ",", j, ") = ", identity(i,j)
                        stop
                    end if
                end do
            end do
            
          
        end subroutine invert_matrix



end module utils