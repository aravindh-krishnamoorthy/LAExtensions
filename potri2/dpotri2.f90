!###############################################################################
! This file is a part of the package: MatrixAlgorithms
! Released under the MIT license, see LICENSE file for details.
! Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
!###############################################################################

SUBROUTINE DPOTRI2(UPLO, N, A, LDA, INFO)

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N
    DOUBLE PRECISION   A( LDA, * )

    INTEGER            ILAENV
    EXTERNAL           ILAENV

    integer            NB
    DOUBLE PRECISION   V(N)
    DOUBLE PRECISION   ONE, ZERO
    PARAMETER ( ONE = 1.0, ZERO = 0.0 )

    ! Reuse ILAENV for DTRITRI
    NB = ILAENV(1, 'DTRITRI', UPLO//'N', N, -1, -1, -1)
    IF (NB.LE.1 .OR. NB.GE.N) THEN

        ! Scalar version
        CALL DPOTRI2S(UPLO, N, A, LDA, INFO)

    ELSE

        ! Future block version
        CALL DPOTRI2S(UPLO, N, A, LDA, INFO)

    END IF

    INFO = 0
    RETURN
END
