module MatrixAlgorithms

using LinearAlgebra

########################################
# potri:
# Positive definite matrix inverse
########################################
include("../potri2/potri2.jl")
########################################
# trsr:
# Matrix square root
########################################
include("../trsr/trsr.jl")
########################################
# ordschur:
# Schur matrix reordering
########################################
include("../ordschur/ordschur.jl")

end
