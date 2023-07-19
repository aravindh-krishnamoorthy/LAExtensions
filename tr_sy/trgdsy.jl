using LinearAlgebra

# Algorithm: 'Matrix Inversion Using Cholesky Decomposition', Aravindh Krishnamoorthy, Deepak Menon, arXiv:1111.4144.
function trgdsy!(uplo::Char, RL::AbstractMatrix{T}, d::AbstractVector{T}) where {T}
    N = size(RL,1)
    for i=1:N
        RL[i,i], d[i] = d[i], RL[i,i]
    end
    if uplo == 'U'
        for j=1:N
            for i=j+1:N
                RL[i,j] = 0
            end
        end
        for j=N:-1:1
            # k = N,...,j+1
            for k=N:-1:j+1
                for i=1:j
                    RL[j,i] = RL[j,i] - RL[i,k]*RL[k,j]
                end
            end
            # k = j,...,1
            for k=j:-1:1
                RL[j,k] = conj(RL[j,k]/d[k])
                for i=1:k-1
                    RL[j,i] = RL[j,i] - RL[i,k]*conj(RL[j,k])
                end
            end
        end
        RL = Hermitian(RL, :L)
    else # if uplo == 'L'
        for j=1:N
            for i=1:j-1
                RL[i,j] = 0
            end
        end
        for j=1:N
            for k=1:j-1
                for i=j:N
                    RL[j,i] = RL[j,i] - RL[i,k]*RL[k,j]
                end
            end
            for k=j:N
                RL[j,k] = conj(RL[j,k]/d[k])
                for i=k+1:N
                    RL[j,i] = RL[j,i] - RL[i,k]*conj(RL[j,k])
                end
            end
        end
        RL = Hermitian(RL, :U)
    end
    return RL
end
