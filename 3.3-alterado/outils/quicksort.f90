SUBROUTINE quick_sort_up(list, n, order)

  ! Sort routine to arrange array elements from smallest to largest.
  ! Grabbed from a millers web site: http://users.bigpond.net.au/amiller

  ! Quick sort routine from:
  ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
  ! Modified by Alan Miller to include an associated integer array which gives
  ! the positions of the elements in the original order.

  implicit none

  integer,                intent(in)    :: n
  real,    dimension (n), intent(inout) :: list
  integer, dimension (n), intent(out)   :: order

  ! LOCAL VARIABLES
  integer :: i

  do i = 1, n
     order(i) = i
  end do

  call quick_sort_up_1(1, n)

CONTAINS

  RECURSIVE SUBROUTINE quick_sort_up_1(left_end, right_end)

    integer, intent(in) :: left_end, right_end

    ! LOCAL VARIABLES
    integer             :: i, j, itemp
    real                :: reference, temp
    integer, parameter  :: max_simple_sort_size = 6

    IF (right_end < left_end + max_simple_sort_size) THEN
       ! Use interchange sort for small lists
       CALL interchange_sort_up(left_end, right_end)

    ELSE
       ! Use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1

       DO
          ! Scan list from left end until element >= reference is found
          DO
             i = i + 1
             IF (list(i) >= reference) EXIT
          END DO
          ! Scan list from right end until element <= reference is found
          DO
             j = j - 1
             IF (list(j) <= reference) EXIT
          END DO


          IF (i < j) THEN
             ! Swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          ELSE IF (i == j) THEN
             i = i + 1
             EXIT
          ELSE
             EXIT
          END IF
       END DO

       IF (left_end < j)  CALL quick_sort_up_1(left_end, j)
       IF (i < right_end) CALL quick_sort_up_1(i, right_end)
    END IF

  END SUBROUTINE quick_sort_up_1


  SUBROUTINE interchange_sort_up(left_end, right_end)

    INTEGER, INTENT(IN) :: left_end, right_end

    ! LOCAL VARIABLES
    INTEGER             :: i, j, itemp
    REAL                :: temp

    DO i = left_end, right_end - 1
       DO j = i+1, right_end
          IF (list(i) > list(j)) THEN
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          END IF
       END DO
    END DO

  END SUBROUTINE interchange_sort_up

END SUBROUTINE quick_sort_up



!===============================================================================



SUBROUTINE quick_sort_down(list, n, order)

  ! Sort routine to arrange array elements from largest to smallest.
  ! Grabbed from a millers web site: http://users.bigpond.net.au/amiller

  ! Quick sort routine from:
  ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
  ! Modified by Alan Miller to include an associated integer array which gives
  ! the positions of the elements in the original order.

  implicit none

  integer,                intent(in)    :: n
  real,    dimension (n), intent(inout) :: list
  integer, dimension (n), intent(out)   :: order

  ! LOCAL VARIABLES
  integer :: i

  do i = 1, n
     order(i) = i
  end do

  call quick_sort_down_1(1, n)

CONTAINS

  RECURSIVE SUBROUTINE quick_sort_down_1(left_end, right_end)

    integer, intent(in) :: left_end, right_end

    ! LOCAL VARIABLES
    integer             :: i, j, itemp
    real                :: reference, temp
    integer, parameter  :: max_simple_sort_size = 6

    IF (right_end < left_end + max_simple_sort_size) THEN
       ! Use interchange sort for small lists
       CALL interchange_sort_down(left_end, right_end)

    ELSE
       ! Use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1

       DO
          ! Scan list from left end until element <= reference is found
          DO
             i = i + 1
             IF (list(i) <= reference) EXIT
          END DO
          ! Scan list from right end until element >= reference is found
          DO
             j = j - 1
             IF (list(j) >= reference) EXIT
          END DO


          IF (i < j) THEN
             ! Swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          ELSE IF (i == j) THEN
             i = i + 1
             EXIT
          ELSE
             EXIT
          END IF
       END DO

       IF (left_end < j)  CALL quick_sort_down_1(left_end, j)
       IF (i < right_end) CALL quick_sort_down_1(i, right_end)
    END IF

  END SUBROUTINE quick_sort_down_1


  SUBROUTINE interchange_sort_down(left_end, right_end)

    INTEGER, INTENT(IN) :: left_end, right_end

    ! LOCAL VARIABLES
    INTEGER             :: i, j, itemp
    REAL                :: temp

    DO i = left_end, right_end - 1
       DO j = i+1, right_end
          IF (list(i) < list(j)) THEN
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          END IF
       END DO
    END DO

  END SUBROUTINE interchange_sort_down

END SUBROUTINE quick_sort_down
