!------------------------------------------------------------------------------
!
! MODULE: mylib
!
!> @author
!> Andrzej Sołtysik
!
! DESCRIPTION: 
!> Matrix multiplication and gaussian elimination, both non-parallelized
!
!------------------------------------------------------------------------------
module mylib
implicit none

contains

!-----
! DESCRIPTION:
!> Multiplies two matrices
!
!> @param[in] a First matrix
!> @param[in] b Second matrix
!> @param[in] a_rows number of rows in first matrix
!> @param[in] a_cols number of columns in first matrix
!> @param[in] b_rows number of rows in second matrix
!> @param[in] b_cols number of columns in second matrix
!> @param[out] ret return matrix
!-----
subroutine matmul(a, b, a_rows, a_cols, b_rows, b_cols, ret)
    INTEGER :: a_rows, a_cols, b_rows, b_cols
    REAL (kind = 4) :: a(a_rows, a_cols)
    REAL (kind = 4) :: b(b_rows, b_cols)
    REAL (kind = 4) :: ret(a_rows, b_cols)

    !f2py intent(in) :: a, b, a_rows, a_cols, b_rows, b_cols
    !f2py intent(out) :: ret

    INTEGER :: r, c, i

    a_rows = size(a, 1)
    a_cols = size(a, 2)
    b_cols = size(b, 2)

    ret = 0

    do r = 1, a_rows
        do c = 1, b_cols
            do i = 1, a_cols
                ret(r, c) = ret(r, c) + a(r, i) * b(i, c)
            enddo
        enddo
    enddo
    
end subroutine matmul

!-----
! DESCRIPTION:
!> Performs gaussian elimination
!
!> @param[inout] a Two-dimensional (n x n) matrix of coefficents
!> @param[inout] x (1 x n) array of right-hand-side values
!> @param[in] n size of first matrix
!-----
subroutine gauss(a, x, n)
    INTEGER (KIND = 8), intent(in) :: n
    REAL (KIND = 8), intent(inout) :: a(n, n), x(n)

    REAL (KIND = 8) :: ratio
    INTEGER (KIND = 8) :: i, j

    do i = 1, n
        x(i) = x(i) / a(i, i)
        a(:, i) = a(:, i) / a(i, i)
        do j = 1, n
            if((i .NE. j) .AND. (ABS(a(i, i) - 0) > 1d-6)) then
                ratio = a(i, j) / a(i, i)
                a(:,j) = a(:, j) - ratio * a(:, i)
                x(j) = x(j) - ratio * x(i)
            endif

        enddo
    enddo
end subroutine gauss

end module mylib