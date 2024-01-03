subroutine dtrgdsy(uplo,n,RL,d)
implicit none
integer n
character uplo
double precision RL(n,n),d(n)
!========================================
integer i,j,k
double precision tmp
do i = 1,n
    tmp = d(i)
    d(i) = RL(i,i)
    RL(i,i) = tmp
end do
if (uplo.EQ.'U') then
    do j = 1,n
        do i = j+1,n
            RL(i,j) = 0
        end do
    end do
    do j = n,1,-1
        do k = n,j+1,-1
            do i = 1,j
                RL(j,i) = RL(j,i) - RL(i,k)*RL(k,j)
            end do
        end do
        do k = j,1,-1
            RL(j,k) = RL(j,k)/d(k)
            do i = 1,k-1
                RL(j,i) = RL(j,i) - RL(i,k)*RL(j,k)
            end do
        end do
    end do
else
    do j = 1,n
        do i= 1,j-1
            RL(i,j) = 0
        end do
    end do
    do j = 1,n
        do k = 1,j-1
            do i = j,n
                RL(j,i) = RL(j,i) - RL(i,k)*RL(k,j)
            end do
        end do
        do k = j,n
            RL(j,k) = RL(j,k)/d(k)
            do i = k+1,n
                RL(j,i) = RL(j,i) - RL(i,k)*RL(j,k)
            end do
        end do
    end do
end if
return
end
