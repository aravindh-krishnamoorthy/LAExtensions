################################################################################
# Invert a positive definite matrix using the algorithm in
# "Matrix Inversion Using Cholesky Decomposition", by Aravindh Krishnamoorthy
#   and Deepak Menon, arXiv:1111.4144.
################################################################################
# Fortran version
################################################################################
function dpotri2!(uplo::Char, X::AbstractMatrix{T}; rl::Bool=true) where {T}
    # DPOTRI2(UPLO, N, A, LDA, INFO)
    N = size(X,1)
    INFO = Ref{Int64}()
    if rl == true
        lib = Libdl.dlopen("./potri2.so")
        dpotri2 = Libdl.dlsym(lib, :dpotri2_)
        ccall(dpotri2, Cvoid,
            (Ref{UInt8}, Ref{Int64}, Ptr{Float64}, Ref{Int64}, Ptr{Int64}, Clong),
            uplo, N, X, N, INFO, 1)
        Libdl.dlclose(lib)
    else
        ccall((:dpotri2_, "./potri2.so"), Cvoid,
            (Ref{UInt8}, Ref{Int64}, Ptr{Float64}, Ref{Int64}, Ptr{Int64}, Clong),
            uplo, N, X, N, INFO, 1)
    end    
    return X
end
