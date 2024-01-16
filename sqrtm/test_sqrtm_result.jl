using LinearAlgebra
using BenchmarkTools

for iter = 1:1000
    for f in [MatrixAlgorithms._sqrt_quasi_triu!, MatrixAlgorithms.d_sqrt_quasi_triu!]
    println("$f $iter/1000...")
    for N in [8, 64, 1000]
            # Real _sqrt_quasi_triu!
            X = randn(N,N)
            X1 = MatrixAlgorithms.sqrtm(X, f)
            @assert X1*X1 ≈ X
            X1 = MatrixAlgorithms.sqrtm(X+X', f)
            @assert X1*X1 ≈ X+X'
            X1 = MatrixAlgorithms.sqrtm(X+N*I(N), f)
            @assert X1*X1 ≈ X+N*I(N)
            X1 = MatrixAlgorithms.sqrtm(X+X'+N*I(N), f)
            @assert X1*X1 ≈ X+X'+N*I(N)
            # Complex _sqrt_quasi_triu!
            X = complex.(randn(N,N),randn(N,N))
            X1 = MatrixAlgorithms.sqrtm(X, f)
            @assert X1*X1 ≈ X
            X1 = MatrixAlgorithms.sqrtm(X+X', f)
            @assert X1*X1 ≈ X+X'
            X1 = MatrixAlgorithms.sqrtm(X+N*I(N), f)
            @assert X1*X1 ≈ X+N*I(N)
            X1 = MatrixAlgorithms.sqrtm(X+X'+N*I(N), f)
            @assert X1*X1 ≈ X+X'+N*I(N)
        end
    end
end
