using LinearAlgebra
using BenchmarkTools

for iter = 1:1000
    println("$iter/1000...")
    for N in [8, 64, 1000]
            # Real trsr!
            X = randn(N,N)
            X1 = MatrixAlgorithms.sqrtm(X)
            @assert X1*X1 ≈ X
            X1 = MatrixAlgorithms.sqrtm(X+X')
            @assert X1*X1 ≈ X+X'
            X1 = MatrixAlgorithms.sqrtm(X+N*I(N))
            @assert X1*X1 ≈ X+N*I(N)
            X1 = MatrixAlgorithms.sqrtm(X+X'+N*I(N))
            @assert X1*X1 ≈ X+X'+N*I(N)
            # Real trsrf!
            X = randn(N,N)
            X1 = MatrixAlgorithms.sqrtm(X, MatrixAlgorithms.trsrf!)
            @assert X1*X1 ≈ X
            X1 = MatrixAlgorithms.sqrtm(X+X', MatrixAlgorithms.trsrf!)
            @assert X1*X1 ≈ X+X'
            X1 = MatrixAlgorithms.sqrtm(X+N*I(N), MatrixAlgorithms.trsrf!)
            @assert X1*X1 ≈ X+N*I(N)
            X1 = MatrixAlgorithms.sqrtm(X+X'+N*I(N), MatrixAlgorithms.trsrf!)
            @assert X1*X1 ≈ X+X'+N*I(N)
            # Complex trsr!
            X = complex.(randn(N,N),randn(N,N))
            X1 = MatrixAlgorithms.sqrtm(X)
            @assert X1*X1 ≈ X
            X1 = MatrixAlgorithms.sqrtm(X+X')
            @assert X1*X1 ≈ X+X'
            X1 = MatrixAlgorithms.sqrtm(X+N*I(N))
            @assert X1*X1 ≈ X+N*I(N)
            X1 = MatrixAlgorithms.sqrtm(X+X'+N*I(N))
            @assert X1*X1 ≈ X+X'+N*I(N)
            # Complex trsrf!
            X = complex.(randn(N,N),randn(N,N))
            X1 = MatrixAlgorithms.sqrtm(X, MatrixAlgorithms.trsrf!)
            @assert X1*X1 ≈ X
            X1 = MatrixAlgorithms.sqrtm(X+X', MatrixAlgorithms.trsrf!)
            @assert X1*X1 ≈ X+X'
            X1 = MatrixAlgorithms.sqrtm(X+N*I(N), MatrixAlgorithms.trsrf!)
            @assert X1*X1 ≈ X+N*I(N)
            X1 = MatrixAlgorithms.sqrtm(X+X'+N*I(N), MatrixAlgorithms.trsrf!)
            @assert X1*X1 ≈ X+X'+N*I(N)
        end
    end
end
