!###############################################################################
! This file is a part of the package: MatrixAlgorithms
! Released under the MIT license, see LICENSE file for details.
! Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
!###############################################################################

SUBROUTINE ZTRSR(N, A, LDA, INFO)

    INTEGER            INFO, LDA, N
    DOUBLE PRECISION   A( LDA, * )

    EXTERNAL            DGEMM, ZDOTU
    COMPLEX*16          ZDOTU

    INTEGER             I
    DOUBLE PRECISION    ONE, ZERO
    PARAMETER ( ONE = 1.0, ZERO = 0.0 )

    INFO = 0
    RETURN
END
