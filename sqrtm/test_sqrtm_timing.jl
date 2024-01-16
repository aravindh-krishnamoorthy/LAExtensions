using LinearAlgebra
using BenchmarkTools

for N in [8, 64, 1000]
    # Real
    println("$N: Real/No structure")
    X = randn(N,N)
    @btime sqrt($X)
    @btime MatrixAlgorithms.sqrtm($X)
    @btime MatrixAlgorithms.sqrtm($X, MatrixAlgorithms.d_sqrt_quasi_triu!)
    println("$N: Real/Symmetric")
    @btime sqrt($X+$X')
    @btime MatrixAlgorithms.sqrtm($X+$X')
    @btime MatrixAlgorithms.sqrtm($X+$X', MatrixAlgorithms.d_sqrt_quasi_triu!)
    println("$N: Real/No structure/PD")
    @btime sqrt($X+$N*I($N))
    @btime MatrixAlgorithms.sqrtm($X+$N*I($N))
    @btime MatrixAlgorithms.sqrtm($X+$N*I($N), MatrixAlgorithms.d_sqrt_quasi_triu!)
    println("$N: Real/Symmetric/PD")
    @btime sqrt($X+$X'+$N*I($N))
    @btime MatrixAlgorithms.sqrtm($X+$X'+$N*I($N))
    @btime MatrixAlgorithms.sqrtm($X+$X'+$N*I($N), MatrixAlgorithms.d_sqrt_quasi_triu!)
    # Complex
    println("$N: Complex/No structure")
    X = complex.(randn(N,N),randn(N,N))
    @btime sqrt($X)
    @btime MatrixAlgorithms.sqrtm($X)
    println("$N: Complex/Hermitian")
    @btime sqrt($X+$X')
    @btime MatrixAlgorithms.sqrtm($X+$X')
    println("$N: Complex/No structure/PD")
    @btime sqrt($X+$N*I($N))
    @btime MatrixAlgorithms.sqrtm($X+$N*I($N))
    println("$N: Complex/Hermitian/PD")
    @btime sqrt($X+$X'+$N*I($N))
    @btime MatrixAlgorithms.sqrtm($X+$X'+$N*I($N))
end
