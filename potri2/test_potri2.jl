using LinearAlgebra
using BenchmarkTools

N = 8
X = rand(N, N)
X = X*X'
U = Matrix(cholesky(X).U)
X1 = MatrixAlgorithms.potri2!(copy(U))
display(norm(X1 - inv(X)))
X1 = MatrixAlgorithms.dpotri2!(copy(U))
display(norm(X1 - inv(X)))
# Timing
@btime LAPACK.potri!('U', copy(U))
@btime MatrixAlgorithms.dpotri2!(copy(U))

X = complex.(rand(N, N), randn(N, N))
X = X*X'
U = Matrix(cholesky(X).U)
X1 = MatrixAlgorithms.potri2!(U)
display(norm(X1 - inv(X)))
