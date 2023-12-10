subroutine dtrtrsy(uplo,n,RL,S)
implicit none
integer n
character uplo
double precision RL(n,n),S(n,n)
!========================================
integer i,j,k
if (uplo.EQ.'U') then
    do concurrent (j = n:1:-1)
        do concurrent (k = n:j+1:-1)
            do concurrent (i = 1:j)
                S(i,j) = S(i,j) - RL(i,k)*S(j,k)
            end do
        end do
        do concurrent (k = j:1:-1)
            S(k,j) = S(k,j)/RL(k,k)
            do concurrent (i = 1:k-1)
                S(i,j) = S(i,j) - RL(i,k)*S(k,j)
            end do
        end do
    end do
else
    do j = 1,n
        do concurrent (k = 1:j-1)
            do concurrent (i = j:n)
                S(i,j) = S(i,j) - RL(i,k)*S(j,k)
            end do
        end do
        do concurrent (k = j:n)
            S(k,j) = S(k,j)/RL(k,k)
            do concurrent (i = k+1:n)
                S(i,j) = S(i,j) - RL(i,k)*S(k,j)
            end do
        end do
    end do
end if
return
end
