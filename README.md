# MatrixAlgorithms

_Staging area for matrix algebra algorithms in **M**ATLAB, **J**ulia, **F**ortran, and **C**/C++_

<div align="center">

  | Directory | Language | Description | Target | Development Stage |
  |---|---|---|---|---|
  | [posi](https://github.com/aravindh-krishnamoorthy/MatrixAlgorithms/tree/main/posol) | M-J-F-C | Efficient inversion of positive definite matrices | LAPACK | Partial |
  | [posp](https://github.com/aravindh-krishnamoorthy/MatrixAlgorithms/tree/main/posol) | M-J-F-C | Efficient solution for `PD_X` in `PD_A*PD_X = PD_B`  | LAPACK | Partial |
  | [sqrtm](https://github.com/aravindh-krishnamoorthy/MatrixAlgorithms/tree/main/sqrtm) | J-F | Efficient square root of matrices |  | Partial |

</div>

---

## Julia Usage
This package can be loaded under Julia as follows:
```julia
julia> include("/<path>/MatrixAlgorithms/src/MatrixAlgorithms.jl")
```
