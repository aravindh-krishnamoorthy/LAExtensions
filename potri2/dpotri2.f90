!###############################################################################
! This file is a part of the package: MatrixAlgorithms
! Released under the MIT license, see LICENSE file for details.
! Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
!###############################################################################

SUBROUTINE DPOTRI2(UPLO, N, A, LDA, INFO)

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N
    DOUBLE PRECISION   A( LDA, * )

    EXTERNAL           DPOTRI2S
    EXTERNAL           DTRTRI, DGEMM

    INTEGER            NB
    DOUBLE PRECISION   ONE, ZERO
    ! Configurable block size
    PARAMETER          (NB = 2)
    DOUBLE PRECISION   V(NB,NB)

    IF (NB.LE.1 .OR. NB.GE.N) THEN
        ! Scalar version
        CALL DPOTRI2S(UPLO, N, A, LDA, INFO)
    ELSE
        DO J = 1, N, NB
            JB = MIN(NB, N-J+1)
            CALL DTRTRI(UPLO, 'N', JB, A(J,J), LDA, INFO)
            DO J1 = J, J+JB-1
                DO I = J1+1, N
                    A(I,J1) = 0
                END DO
            END DO
        END DO
        NN = ((N-1)/NB)*NB + 1
        DO J = NN, 1, -NB
            JB = MIN(NB, N-J+1)
            ! [FIRST PART HERE]
            DO K = J, 1, -JB
                CALL DGEMM('N', 'T', JB, JB, JB, 1.0D0, A(J,K), LDA, A(K,K), LDA, 0.0D0, V, NB)
                A(J:J+JB-1, K:K+JB-1) = V(1:JB, 1:JB)
            END DO
        END DO
    END IF
    INFO = 0
    RETURN
END
