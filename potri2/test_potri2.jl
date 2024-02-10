################################################################################
# This file is a part of the package: MatrixAlgorithms
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
################################################################################

# using MKL
using LinearAlgebra
using BenchmarkTools

for N in [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
    println("N=$N...")

    # Real
    X = rand(N, N)
    X = X*X'
    U = Matrix(cholesky(X).U)
    L = Matrix(cholesky(X).L)
    X0 = LAPACK.potri!('U', copy(U))
    X1 = MatrixAlgorithms.potri2!('U', copy(U))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.potri2!('L', copy(L))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.dpotri2!('U', copy(U))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.dpotri2!('L', copy(L))
    display(norm(triu(X1) - X0))
    # Timing
    @btime LAPACK.potri!('U', copy($U)) ;
    @btime MatrixAlgorithms.potri2!('U', copy($U)) ;
    @btime MatrixAlgorithms.potri2!('L', copy($L)) ;
    @btime MatrixAlgorithms.dpotri2!('U', copy($U); rl=false) ;
    @btime MatrixAlgorithms.dpotri2!('L', copy($L); rl=false) ;
end
