using LinearAlgebra

# Algorithm: 'Matrix Inversion Using Cholesky Decomposition', Aravindh Krishnamoorthy, Deepak Menon, arXiv:1111.4144.
function trtrsy!(uplo::Char, RL::AbstractMatrix{T}, S::AbstractMatrix{T}) where {T}
    N = size(RL,1)
    if uplo == 'U'
        for j=N:-1:1
            # k = N,...,j+1
            for k=N:-1:j+1
                for i=1:j
                    S[i,j] = S[i,j] - RL[i,k]*conj(S[j,k])
                end
            end
            # k = j,...,1
            for k=j:-1:1
                S[k,j] = S[k,j]/RL[k,k]
                for i=1:k-1
                    S[i,j] = S[i,j] - RL[i,k]*S[k,j]
                end
            end
        end
        S = Hermitian(S, :U)
    else # if uplo == 'L'
        for j=1:N
            # k = i+1,...,N
            for k=1:j-1
                for i=j:N
                    S[i,j] = S[i,j] - RL[i,k]*conj(S[j,k])
                end
            end
            # k = 1,...,i
            for k=j:N
                S[k,j] = S[k,j]/RL[k,k]
                for i=k+1:N
                    S[i,j] = S[i,j] - RL[i,k]*S[k,j]
                end
            end
        end
        S = Hermitian(S, :L)
    end
    return S
end

