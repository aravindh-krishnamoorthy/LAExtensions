using LinearAlgebra

N = 8
X = rand(N, N)
X = X*X'
U = Matrix(cholesky(X).U)
X1 = MatrixAlgorithms.potri2!(U)
display(norm(X1 - inv(X)))

X = complex.(rand(N, N), randn(N, N))
X = X*X'
U = Matrix(cholesky(X).U)
X1 = MatrixAlgorithms.potri2!(U)
display(norm(X1 - inv(X)))
