using Random
using LinearAlgebra
include("../src/MatrixAlgorithms.jl")

rng = MersenneTwister(555);
A = randn(rng,64,64)
A = A*A'

T = Matrix(cholesky(A).U) ;

T1 = LAPACK.potri!('U', copy(T))
T1 = T1+triu(T1,1)'
display(T1)
T2 = MatrixAlgorithms.dpotri2!('U', copy(T))
display(T2)
display(T1-T2)
norm(T1-T2)
