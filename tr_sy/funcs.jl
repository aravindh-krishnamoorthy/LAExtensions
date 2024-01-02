using LinearAlgebra

# Include lower level files
include("trgdsy.jl")

# Select default implementation
inv(A::Union{Hermitian{T,S}, Symmetric{T,S}}) where {T,S} = inv(A, ReferenceImpl())

# Reference Implementation
function inv(A::Union{Hermitian{T,S}, Symmetric{T,S}}, ::ReferenceImpl) where {T,S}
    LU = lu(A, NoPivot())
    return trgdsy!('U', LU.U, ones(size(A,1)))
end

# Select default implementation
inv(C::CholeskyPivoted{T}) where {T} = inv(C, ReferenceImpl())

# Reference Implementation
function inv(C::CholeskyPivoted{T}, ::ReferenceImpl) where {T}
    R = copy(C.U)
    p = invperm(C.p)
    return trgdsy!('U', Matrix{T}(R), 1 ./diag(R))[p,p]
end
