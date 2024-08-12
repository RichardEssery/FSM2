!-----------------------------------------------------------------------
! Solve matrix equation Ax = b for x by LU decomposition
!-----------------------------------------------------------------------
subroutine LUDCMP(N,A,b,x)

implicit none

integer, intent(in) :: &
  N                   ! Number of equations to solve

real, intent(in) :: &
 A(N,N),            & ! Matrix
 b(N)                 ! RHS of matrix equation

real, intent(out) :: &
 x(N)                 ! Solution of matrix equation

integer :: i,ii,imax,j,k,ll,indx(N)

real :: Acp(N,N),aamax,dum,sum,vv(N)

Acp(:,:) = A(:,:)
x(:) = b(:)

do i = 1, N
  aamax = 0
  do j = 1, N
    if (abs(Acp(i,j)) > aamax) aamax = abs(Acp(i,j))
  end do
  if (aamax == 0) return
  vv(i) = 1/aamax
end do

do j = 1, N
  do i = 1, j - 1
    sum = Acp(i,j)
    if (i > 1) then
      do k = 1, i - 1
        sum = sum - Acp(i,k)*Acp(k,j)
      end do
      Acp(i,j) = sum
    end if
  end do
  aamax = 0
  do i = j, N
    sum = Acp(i,j)
    do k = 1, j - 1
      sum = sum - Acp(i,k)*Acp(k,j)
    end do
    Acp(i,j) = sum
    dum = vv(i)*abs(sum)
    if (dum >= aamax) then
      imax = i
      aamax = dum
    end if
  end do
  if (j /= imax)then
    do k = 1, N
      dum = Acp(imax,k)
      Acp(imax,k) = Acp(j,k)
      Acp(j,k) = dum
    end do
    vv(imax) = vv(j)
  end if
  indx(j) = imax
  if (Acp(j,j) == 0) Acp(j,j) = 1e-20
  if (j /= N) then
    dum = 1/Acp(j,j)
    do i = j + 1, N
      Acp(i,j) = Acp(i,j)*dum
    end do
  end if
end do

ii = 0
do i = 1, N
  ll = indx(i)
  sum = x(ll)
  x(ll) = x(i)
  if (ii /= 0)then
    do j = ii, i - 1
      sum = sum - Acp(i,j)*x(j)
    end do
  else if (sum /= 0) then
    ii=i
  end if
  x(i) = sum
end do

do i = N, 1, -1
  sum = x(i)
  do j = i + 1, N
    sum = sum - Acp(i,j)*x(j)
  end do
  x(i) = sum/Acp(i,i)
end do
    
end subroutine LUDCMP
