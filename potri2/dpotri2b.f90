!###############################################################################
! This file is a part of the package: MatrixAlgorithms
! Released under the MIT license, see LICENSE file for details.
! Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
!###############################################################################

SUBROUTINE DPOTRI2B(UPLO, N, A, LDA, INFO)

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N
    DOUBLE PRECISION   A( LDA, * )

    EXTERNAL           DPOTRI2S
    EXTERNAL           DTRTRI, DGEMM

    INTEGER            NB
    ! Configurable block size
    PARAMETER          (NB = 2)
    DOUBLE PRECISION   V(NB,N)
    DOUBLE PRECISION   W(NB,NB)

    DO J = 1, N, NB
        JB = MIN(NB, N-J+1)
        CALL DTRTRI(UPLO, 'N', JB, A(J,J), LDA, INFO)
        V(1:JB,J:J+JB-1) = A(J:J+JB-1,J:J+JB-1)
        DO J1 = J, J+JB-1
            DO I = J1+1, N
                A(I,J1) = 0
            END DO
        END DO
    END DO
    NN = ((N-1)/NB)*NB + 1
    DO J = NN, 1, -NB
        JB = MIN(NB, N-J+1)
        DO I = 1, J, JB
            DO K = NN, J+1, -JB
                CALL DGEMM('T', 'T', JB, JB, JB, -1.0D0, A(K,J), LDA, A(I,K), LDA, 1.0D0, A(J,I), LDA)
            END DO
        END DO
        DO K = J, 1, -JB
            CALL DGEMM('N', 'T', JB, JB, JB, 1.0D0, A(J,K), LDA, V(1,K), NB, 0.0D0, W, NB)
            A(J:J+JB-1, K:K+JB-1) = W(1:JB, 1:JB)
            DO I = 1, K-1, JB
                CALL DGEMM('N', 'T', JB, JB, JB, -1.0D0, A(J,K), LDA, A(I,K), LDA, 1.0D0, A(J,I), LDA)
            END DO
        END DO
    END DO
    DO I = 1, N
        DO J = I+1, N
            A(I,J) = A(J,I)
        END DO
    END DO
    INFO = 0
    RETURN
END
