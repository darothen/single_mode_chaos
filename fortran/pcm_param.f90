!-----------------
module param

    implicit none
    integer, parameter :: dp = kind(1.d0)
    character(*), parameter :: pcm_table = "SM_OLS_4.ascii"

    ! SET THESE PARAMETERS
    logical            :: LOG_SMAX = .true.
    integer, parameter :: MAX_ORDER = 5
    !integer, parameter :: NVARS = 8
    integer, parameter :: max_vars = 40
    integer            :: NTERMS ! to be set when reading the coeffs

    integer                               :: NVARS
    character(len=10), dimension(max_vars) :: var_names
    real (kind=dp), dimension(max_vars)   :: lower_bnds, upper_bnds
    logical, dimension(max_vars)          :: log_var
    logical                               :: PRNT_DBG 
    namelist /pcm/ NVARS, var_names, log_var, lower_bnds, upper_bnds, PRNT_DBG

    real (kind=dp), dimension(:), allocatable :: coeffs
    integer, dimension(:,:), allocatable :: orders 
    
contains

    subroutine maxsat( xin, Smax )
        real (kind=dp), dimension(:) :: xin
        real (kind=dp), intent(out) :: Smax

        real (kind=dp), dimension(NVARS, 0:MAX_ORDER) :: ortho_poly_evals
        real (kind=dp), dimension(NVARS) :: zin
        integer :: i, j

        ! 1) For all the log variables, apply the correct
        ! transformation (might need to be tweaked by hand)
        do i = 1, NVARS
            if (log_var(i)) then
                xin(i) = log10(xin(i))
            end if
        end do

        ! 2) Project into z-space
        call project ( xin, zin )
        if (PRNT_DBG) then
            do i = 1, NVARS
                print *, xin(i), zin(i)
            end do
        end if

        ! 3) Evaluate necessary orthogonal polynomials
        call p_polynomial_value ( NVARS, MAX_ORDER, zin, ortho_poly_evals )
        if (PRNT_DBG) then
        print *, "ORTHO_POLY_EVALS"
            do i = 1, NVARS
                do j = 0, MAX_ORDER
                    print *, i, j, ortho_poly_evals(i, j)
                end do
            end do
        end if

        ! 4) Evaluate the chaos expansion using the stored 
        ! orthogonal polynomial evaluations
        call poly_eval( orders, coeffs, ortho_poly_evals, Smax )

        ! 5) Post-process to the value we care about
        if ( LOG_SMAX ) then
            Smax = 10.d0**Smax
        end if

        return

    end subroutine maxsat

    subroutine poly_finish

        deallocate( coeffs, orders )

    end subroutine poly_finish

    subroutine poly_init

        integer :: iostat, i, j, counter
        real (kind=dp) :: row_coeff
        real (kind=dp), dimension(NVARS) :: row_orders

        open(unit=10, file=pcm_table, status='OLD', action='READ', &
             iostat=iostat)
        if (iostat /= 0) then
            print *, "There was a problem opening" // pcm_table // "...terminating"
            stop
        end if 

        ! First, skim over the saved PCM coefficients/orders 
        counter = 0
        inspect_file : do
            read (10, *, iostat=iostat) 
            if (iostat > 0) then
                print *, "Error reading " // pcm_table
                stop
            else if (iostat < 0) then
                print *, "Finished reading " // pcm_table
                exit
            else ! successfully read a line
                counter = counter + 1
                !print *, counter, row_orders
            end if
        end do inspect_file

        rewind (10)
        NTERMS = counter

        ! Appropriate enough memory to store the table
        allocate( coeffs(NTERMS), orders(nterms,NVARS) )

        ! Now we can actually read and save it
        print *, "actual read loop - nterms= ", NTERMS
        read_file : do i = 1, NTERMS
            read (10, *, iostat=iostat) row_coeff, row_orders
            if (iostat > 0) then
                print *, "Error reading " // pcm_table
                stop
            else if (iostat < 0) then
                print *, "Finished reading " // pcm_table
                exit
            else ! successfully read a line
                coeffs(i) = row_coeff
                do j = 1, NVARS
                    orders(i,j) = int(row_orders(j))
                end do
            end if
        end do read_file

        if (PRNT_DBG) then
            print *, "POLY_INIT RESULTS"
            print *, "COEFFS"
            do i = 1, NTERMS
                print *, i, coeffs(i)
            end do
            print *, "ORDERS"
            do i = 1, NTERMS
                do j = 1, NVARS
                    print *, i, j, orders(i, j)
                end do
            end do
            print *, "----------"
        end if

        close (unit=10)

    end subroutine

    subroutine poly_eval( orders, coeffs, zin, eval )
        integer, dimension(NTERMS,NVARS), intent(in) :: orders
        real (kind=dp), dimension(NVARS, 0:MAX_ORDER), intent(in) :: zin
        real (kind=dp), dimension(NVARS), intent(in) :: coeffs
        real (kind=dp), intent(out) :: eval

        real (kind=dp) :: row_eval
        integer :: row, col

        if (PRNT_DBG) then
            print *, "poly_eval sr", NTERMS, NVARS
            print *, "orders", size(orders, dim=1), size(orders, dim=2)
            print *, "zin", size(zin, dim=1), size(zin, dim=2) 
        end if

        eval = 0d0
        do row = 1, NTERMS
            row_eval = 1d0
            do col = 1, NVARS
                if (PRNT_DBG) then
                    print *, row, col, orders(row, col), zin(col, orders(row, col))
                end if
                row_eval = row_eval*zin(col, orders(row, col))
            end do
            if (PRNT_DBG) then
                print *, row, row_eval, coeffs(row)
            end if
            eval = eval + row_eval*coeffs(row)
        end do

        return 

    end subroutine poly_eval

    subroutine project ( xin, zout )
        real (kind=dp), dimension(NVARS), intent(in) :: xin
        real (kind=dp), dimension(NVARS), intent(out) :: zout

        real (kind=dp) :: x, z
        integer :: i

        do i = 1, NVARS
            x = xin(i)
            z = uni_to_uni(x, lower_bnds(i), upper_bnds(i))
            zout(i) = z
        end do

        return

    end subroutine

    real (kind=dp) function uni_to_uni (x, ai, bi)
        real (kind=dp), intent(in) :: x, ai, bi
        ! Hard code af, bf to [-1, 1]
        real (kind=dp) :: af=-1d0, bf=1d0

        if (x < ai) then
            uni_to_uni = af
        else if (x > bi) then
            uni_to_uni = bf
        else
            uni_to_uni = ((bf-af)/(bi-ai))*(x-ai) + af
        end if

    end function uni_to_uni

    subroutine p_polynomial_value ( m, n, x, v )

    !*****************************************************************************80
    !
    !! P_POLYNOMIAL_VALUE evaluates the Legendre polynomials P(n,x).
    !
    !  Discussion:
    !
    !    P(n,1) = 1.
    !    P(n,-1) = (-1)^N.
    !    | P(n,x) | <= 1 in [-1,1].
    !
    !    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
    !    quadrature of the integral of a function F(X) with weight function 1
    !    over the interval [-1,1].
    !
    !    The Legendre polynomials are orthogonal under the inner product defined
    !    as integration from -1 to 1:
    !
    !      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX 
    !        = 0 if I =/= J
    !        = 2 / ( 2*I+1 ) if I = J.
    !
    !    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
    !
    !    A function F(X) defined on [-1,1] may be approximated by the series
    !      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
    !    where
    !      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,x) dx.
    !
    !    The formula is:
    !
    !      P(n,x) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
    !
    !  Differential equation:
    !
    !    (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0
    !
    !  First terms:
    !
    !    P( 0,x) =      1
    !    P( 1,x) =      1 X
    !    P( 2,x) = (    3 X^2 -       1)/2
    !    P( 3,x) = (    5 X^3 -     3 X)/2
    !    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
    !    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
    !    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
    !    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
    !    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
    !    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
    !    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
    !
    !  Recursion:
    !
    !    P(0,x) = 1
    !    P(1,x) = x
    !    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n
    !
    !    P'(0,x) = 0
    !    P'(1,x) = 1
    !    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    10 March 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    National Bureau of Standards, 1964,
    !    ISBN: 0-486-61272-4,
    !    LC: QA47.A34.
    !
    !    Daniel Zwillinger, editor,
    !    CRC Standard Mathematical Tables and Formulae,
    !    30th Edition,
    !    CRC Press, 1996.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, the number of evaluation points.
    !
    !    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
    !    Note that polynomials 0 through N will be evaluated.
    !
    !    Input, real ( kind = 8 ) X(M), the evaluation points.
    !
    !    Output, real ( kind = 8 ) V(M,0:N), the values of the Legendre polynomials 
    !    of order 0 through N at the points X.
    !

      integer ( kind = 4 ) m
      integer ( kind = 4 ) n

      integer ( kind = 4 ) i
      real ( kind = dp ) v(m,0:n)
      real ( kind = dp ) x(m)

      if ( n < 0 ) then
        return
      end if

      v(1:m,0) = 1.0D+00

      if ( n < 1 ) then
        return
      end if

      v(1:m,1) = x(1:m)
     
      do i = 2, n
     
        v(1:m,i) = ( real ( 2 * i - 1, kind = 8 ) * x(1:m) * v(1:m,i-1)   &
                   - real (     i - 1, kind = 8 ) *          v(1:m,i-2) ) &
                   / real (     i,     kind = 8 )
     
      end do
     
      return
    
    end subroutine p_polynomial_value

end module param

program main

    use param

    implicit none

    integer :: iostat, i
    real (kind=dp) :: Smax
    real (kind=dp), dimension(:), allocatable :: x
    real (kind=dp) :: n = 850.0
    real (kind=dp) :: mu = 0.05

    open(10,file='pcm.nl', iostat=iostat)
    if (iostat /= 0) then
        print *, "Could not open open 'pcm.nl'"
        stop
    end if 
    read(10, nml=pcm)
    do i = 1, NVARS
        print *, i, var_names(i), log_var(i), lower_bnds(i), upper_bnds(i)
    end do

    allocate( x(NVARS) )

    !x = (/ 1., 0.004, 2.0, 0.507, 0.2, 293., 85000., 0.1 /)
    x = (/ 0., 0., &
           2.0, 0.507, 0.2, 266.87, 91726., 0.1 /)
    x(1) = n
    x(2) = mu

    print *, "Initializing saved PCM..."

    call poly_init

    call maxsat ( x, Smax )

    call poly_finish

    print *, Smax

    !print *, nvars
    !print *, var_names
    !print *, log_var

    deallocate( x )

end program main

