!###############################################################################
! This file is a part of the package: MatrixAlgorithms
! Released under the MIT license, see LICENSE file for details.
! Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
!###############################################################################

SUBROUTINE DPOTRI2B(UPLO, N, A, LDA, INFO)
    IMPLICIT           NONE

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N
    DOUBLE PRECISION   A( LDA, * )

    INTEGER            NB

    ! Configurable block size
    PARAMETER          (NB = 32)

    INFO = 0
    RETURN
END
