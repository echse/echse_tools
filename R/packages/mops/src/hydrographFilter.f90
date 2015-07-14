module sub
contains

subroutine heapsort(arr)
implicit none
!intent(out)
real(8),dimension(:),intent(inout):: arr
!local
real(8):: tmp
integer(4):: i,n
!code
n=size(arr)
do i=(n/2),1,-1
  CALL sift_down(i,n)
end do
do i=n,2,-1
  tmp= arr(1)
  arr(1)= arr(i)
  arr(i)= tmp
  CALL sift_down(1,i-1)
end do
contains
  subroutine sift_down(l,r)
  implicit none
  integer(4),intent(in):: l,r
  integer(4):: j,jold
  real(8):: a
  a=arr(l)
  jold=l
  j=l+l
  do
    if (j .gt. r) then
      exit
    end if
    if (j .lt. r) then
      if (arr(j) .lt. arr(j+1)) then
        j=j+1
      end if
    end if
    if (a .ge. arr(j)) then
      exit
    end if
    arr(jold)=arr(j)
    jold=j
    j=j+j
  end do
  arr(jold)=a
  end subroutine sift_down
end subroutine heapsort


function quantile(array,alpha) result (res)
implicit none
!intent(in)
real(8),dimension(:),intent(in):: array
real(8),intent(in):: alpha
!local
real(8):: res
real(8),dimension(size(array)):: array_sorted
real(8):: pos,posfrac,lower,upper
real(8),parameter:: ZERO= 0.0_8, ONE= 1.0_8
!code
if (alpha .le. ZERO) then
  res= minval(array)
else if (alpha .ge. ONE) then
  res= maxval(array)
else
  array_sorted(:)= array(:)
  call heapsort(array_sorted(:))
  if (size(array_sorted) .eq. 1) then
    res= array_sorted(1)
  else
    pos = ONE + alpha * dble(size(array_sorted) - 1) !Alpha=0 --> pos=1, Alpha=1 --> pos=n
    posfrac= pos - floor(pos)                        !Fractional part of pos
    lower = array_sorted(floor(pos))
    upper = array_sorted(floor(pos)+1)
    res= lower + posfrac * (upper - lower)
  end if
end if
end function quantile

end module sub

!-------------------------------------------------------------------------------

subroutine selectData(vect, vect_len, vect_na, win_len, win_lag, alpha, margins_ok, selected)
  use sub
  implicit none
  ! Intent(in)
  ! Vector with data
  real(8),dimension(vect_len),intent(in):: vect
  integer,intent(in):: vect_len
  real(8),intent(in):: vect_na
  ! Window specification
  integer,intent(in):: win_len
  integer,intent(in):: win_lag
  ! Selection threshold
  real(8),intent(in):: alpha
  ! Is selectio at window margins ok?
  logical,intent(in):: margins_ok
  ! Intent(out)
  logical,dimension(vect_len),intent(out):: selected
  ! Local
  integer:: mini, maxi, i
  logical,dimension(win_len):: local
  ! Code
  ! Initialize result
  selected= .false.
  ! Check errors
  if (win_lag .ge. 1) then
    ! Initialize indices
    mini= 1
    maxi= win_len
    ! Move through vector
    do while (maxi .le. vect_len)
      ! Initialize local window
      local= .false.
      ! Select data
      if (.not. any(vect(mini:maxi) .eq. vect_na)) then
        where(vect(mini:maxi) .gt. quantile(vect(mini:maxi), alpha)) local= .true.
      end if
      ! Undo selection at margins if requested
      if (.not. margins_ok) then
        do i=1,win_len
          if (.not. local(i)) exit
          local(i)= .false.
        end do
        do i=win_len,1,-1
          if (.not. local(i)) exit
          local(i)= .false.
        end do      
      end if
      ! Save result for local window (without unselecting formerly selected data)
      if (any(local)) where (local) selected(mini:maxi)= .true.
      ! Update indices
      mini= mini + win_lag
      maxi= maxi + win_lag
    end do
  end if
end subroutine selectData

