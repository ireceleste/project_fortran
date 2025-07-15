module utils
    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
    real(dp), parameter :: e = exp(1.0_dp)
    
    interface xgetrf
        module procedure dgetrf_interface
    end interface

    interface xgetri
        module procedure dgetri_interface
    end interface

    contains    

        ! Function to convert a string to lowercase
        pure function to_lower(str) result(lower)
            use iso_fortran_env, only: output_unit

            character(len=*), intent(in) :: str
            character(len=len(str)) :: lower
            integer :: i, c
            

            do i = 1, len(str)
                c = ichar(str(i:i))
                if (c >= ichar('A') .and. c <= ichar('Z')) then
                    lower(i:i) = char(c + 32)
                else
                    lower(i:i) = str(i:i)
                end if
            end do

        end function to_lower

        ! Function to read 1,2,3 columns of numbers from a file
        subroutine read_from_file(filename, n, xvals, yvals, yerrs)

            character(len=*), intent(in) :: filename
            integer, intent(out) :: n
            real(dp), allocatable, intent(out) :: xvals(:)
            real(dp), allocatable, intent(out), optional :: yvals(:), yerrs(:)

            integer :: i, ios, iu = 11, line_count
            character(len=256) :: msg
            real(dp) :: dummy
            logical :: exists

            ! Open and count lines
            inquire(file=filename, exist=exists)
            if (.not. exists) then
                write(*, '(/,A,A)') "Error: File '", trim(filename), "' does not exist."
                stop
            end if

            open(newunit=iu, file=filename, status='old', action='read', iostat=ios, iomsg=msg)
            if (ios /= 0) then
                write(*, '(/,A,A)') "Error opening file: ", trim(msg)
                stop
            end if

            line_count = 0
            do
                read(iu, *, iostat=ios) dummy
                if (ios /= 0) exit
                line_count = line_count + 1
            end do
            n = line_count

            allocate(xvals(n))
            if(present(yvals)) then
                allocate(yvals(n))
            end if
            if(present(yerrs)) then
                allocate(yerrs(n))
            end if

            rewind(iu)
            do i = 1, n
                if(present(yvals) .and. present(yerrs)) then
                    read(iu, *, iostat=ios) xvals(i), yvals(i), yerrs(i)
                else if(present(yvals)) then
                    read(iu, *, iostat=ios) xvals(i), yvals(i)
                else
                    read(iu, *, iostat=ios) xvals(i)
                end if




                if (ios /= 0) then
                    write(*, '("Read error at line ", I0)') i
                    stop
                end if
            end do

            close(iu)
            write(*, '("Read ",I0, " data points from file " ,A)')  n, trim(filename)

        end subroutine read_from_file


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
            integer :: n, info, lwork, i, j
            integer, allocatable :: ipiv(:)
            real(dp), allocatable :: work(:), identity(:,:)

            if (size(A, 1) /= size(A, 2)) then
                print *, "Error: invert_matrix expects a square matrix."
                stop
            end if

            n = size(A, 1)
            Ainv = A
            allocate(ipiv(n))

            ! LU decomposition
            call xgetrf(n, n, Ainv, n, ipiv, info) 
            if (info /= 0) then
                write(*,'(A)') "################### Error in LU decomposition ###################"
                write(*, '( " LU decomposition failed. INFO = " , I0)') info
                if(info < 0) then
                   write(*, '( " The ", I0 , "th argument had an illegal value.")') -info
                else if (info > 0) then
                    write(*, '( "The factor U is exactly singular; the ", I0 , "th diagonal element is zero.")') info
                end if
                write(*,'(A)') "#################################################################"
                stop
            end if

            ! Query for optimal workspace size
            lwork = -1 
            allocate(work(1))
            call xgetri(n, Ainv, n, ipiv, work, lwork, info)
            lwork = int(work(1))  
            deallocate(work)
            allocate(work(lwork))

            ! Perform matrix inversion
            call xgetri(n, Ainv, n, ipiv, work, lwork, info)
            if (info /= 0) then
                print *, "################### Error in matrix inversion ###################"
                write(*, '( " Matrix inversion failed. INFO = " , I0)') info
                if(info < 0) then
                   write(*, '( " The ", I0 , "th argument had an illegal value.")') -info
                else if (info > 0) then
                    write(*, '( "The factor U is exactly singular; the ", I0 , "th diagonal element is zero.")') info
                end if
                write(*,'(A)') "#################################################################"
                stop
            end if

            deallocate(ipiv, work)

            ! Verify inversion

            allocate(identity(n,n))
            identity = matmul(A, Ainv)

            do i = 1, n
                do j = 1, n
                    if ( abs(identity(i,j) - merge(1.0_dp, 0.0_dp, i == j)) > 1.0e-10_dp ) then
                        write(*,*) "################### Inversion check failed ###################"
                        write(*, '(" Identity matrix element (",I0,",",I0,") = ",F0.10)')  i,  j, identity(i,j)
                        call exit(1)
                    end if
                end do
            end do
            
            deallocate(identity)
          
        end subroutine invert_matrix



end module utils