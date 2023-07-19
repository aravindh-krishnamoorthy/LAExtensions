module LAPACKExtensions

using LinearAlgebra

const ExternalImpl = Val{1}
const ReferenceImpl = Val{1}

# tr_sy
include("../tr_sy/funcs.jl")

end
