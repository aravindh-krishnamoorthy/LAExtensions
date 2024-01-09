################################################################################
# This file is a part of the package: MatrixAlgorithms
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
################################################################################

using LinearAlgebra
using BenchmarkTools

N = 1000
# Real
X = rand(N, N)
X = X*X'
U = Matrix(cholesky(X).U)
X0 = LAPACK.potri!('U', copy(U))
X1 = MatrixAlgorithms.potri2!('U', copy(U))
display(norm(triu(X1) - X0))
X1 = MatrixAlgorithms.dpotri2!('U', copy(U))
display(norm(triu(X1) - X0))
X1 = MatrixAlgorithms.dpotri2!('L', copy(U'))
display(norm(triu(X1) - X0))
# Timing
@btime LAPACK.potri!('U', copy(U)) ;
@btime MatrixAlgorithms.potri2!('U', copy(U)) ;
@btime MatrixAlgorithms.dpotri2!('U', copy(U)) ;
@btime MatrixAlgorithms.dpotri2!('L', copy(U')) ;

# # Complex
# X = complex.(rand(N, N), randn(N, N))
# X = X*X'
# U = Matrix(cholesky(X).U)
# X0 = LAPACK.potri!('U', copy(U))
# X1 = MatrixAlgorithms.potri2!(copy(U))
# display(norm(triu(X1) - X0))

