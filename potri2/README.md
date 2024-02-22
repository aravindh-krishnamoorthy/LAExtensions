# potri2
An efficient version of LAPACK's `potri` using the algorithm in "Matrix Inversion Using Cholesky Decomposition," by Aravindh Krishnamoorthy and Deepak Menon, [arXiv:1111.4144](https://arxiv.org/abs/1111.4144).

**Note:** This function is only partially implemented.

- For $N \leq 32,$ the scalar version in `dpotri2s.f90` is used, which is partially optimized.
- For $N>32,$ the block version in `dpotri2b.f90` is used, which is not yet optimized.
- Only the IEEE double precision version is currently being implemented. Once complete, the single precision and complex-valued variants will be implemented.

## IEEE double precision

Linux (x86_64-linux-gnu) on 8×11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz for commit [2d93904](https://github.com/aravindh-krishnamoorthy/MatrixAlgorithms/commit/2d93904edbdbd0ce4657899aff0e3d7eb7df8e62)

| N | MKL/U (ns) | MKL/L (ns) | Fortran/U (ns) | Fortran/L (ns) |
| :--- | :--- | :--- | :--- | :--- |
| 2 | 2.7e+02 | 2.6e+02 | 45 | 46 |
| 4 | 4.4e+02 | 4.2e+02 | 73 | 81 |
| 8 | 8.2e+02 | 7.8e+02 | 2.1e+02 | 2.4e+02 |
| 16 | 2.1e+03 | 2e+03 | 1.2e+03 | 1.2e+03 |
| 32 | 4.8e+03 | 4.9e+03 | 5.9e+03 | 5.4e+03 |
| 64 | 1.6e+04 | 2.2e+04 | 4.6e+04 | 3.4e+04 |
| 128 | 9.2e+04 | 1.2e+05 | 4.8e+05 | 2.6e+05 |
| 256 | 5.8e+05 | 6.3e+05 | 4.9e+06 | 2.9e+06 |
