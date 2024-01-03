subroutine dtrtrsy(uplo,n,RL,S)
implicit none
integer n
character uplo
double precision RL(n,n),S(n,n)
!========================================
integer i,j,k
if (uplo.EQ.'U') then
    do j = n,1,-1
        do k = n,j+1,-1
            do i = 1,j
                S(i,j) = S(i,j) - RL(i,k)*S(j,k)
            end do
        end do
        do k = j,1,-1
            S(k,j) = S(k,j)/RL(k,k)
            do i = 1,k-1
                S(i,j) = S(i,j) - RL(i,k)*S(k,j)
            end do
        end do
    end do
else
    do j = 1,n
        do k = 1,j-1
            do i = j,n
                S(i,j) = S(i,j) - RL(i,k)*S(j,k)
            end do
        end do
        do k = j,n
            S(k,j) = S(k,j)/RL(k,k)
            do i = k+1,n
                S(i,j) = S(i,j) - RL(i,k)*S(k,j)
            end do
        end do
    end do
end if
return
end
