#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

int main(int argc, char *argv[]) {
    // Parametri matrica i konstante
    integer n_kvadrat = 4, n = 2, ione = 1; 
    doublereal scalar_one = 1, scalar_mone = -1, scalar_zero = 0;

    // Alokacija memorije
    doublereal *A = malloc(n_kvadrat * n * sizeof(doublereal));
    doublereal *B = malloc(n_kvadrat * n * sizeof(doublereal));
    doublereal *C = malloc(n * n * sizeof(doublereal));
    doublereal *Q, *matrica_Singularne, *U, *Vt, *WORK, *D;

    // Inicijalizacija matrica A i B
    A[0] = 1.2; A[1] = 2.9; A[2] = 5.2; A[3] = 6.8;
    A[4] = 2.1; A[5] = 4.3; A[6] = 6.1; A[7] = 8.1;

    B[0] = 1.0; B[1] = 3.0; B[2] = 5.0; B[3] = 7.0;
    B[4] = 2.0; B[5] = 4.0; B[6] = 6.0; B[7] = 8.0;

    // Množenje matrica A^T i B
    dgemm_("T", "N", &n, &n, &n_kvadrat, 
           &scalar_one, A, &n_kvadrat, B, &n_kvadrat, 
           &scalar_zero, C, &n);

    // Singularna dekompozicija matrice C
    integer LWORK = 5 * n, info;
    matrica_Singularne = malloc(n * sizeof(doublereal));
    U = malloc(n * n * sizeof(doublereal));
    Vt = malloc(n * n * sizeof(doublereal));
    WORK = malloc(LWORK * sizeof(doublereal));

    dgesvd_("A", "A", &n, &n, C, &n, 
            matrica_Singularne, U, &n, Vt, &n, 
            WORK, &LWORK, &info);

    // Računanje matrice Q = U * V^T
    Q = malloc(n * n * sizeof(doublereal));
    dgemm_("N", "N", &n, &n, &n, 
           &scalar_one, U, &n, Vt, &n, 
           &scalar_zero, Q, &n);

    // Ispis matrice Q
    printf("Matrica Q:\n");
    for (integer i = 0; i < n; i++) {
        for (integer j = 0; j < n; j++) {
            printf("%.5lf ", Q[i + n * j]);
        }
        printf("\n");
    }

    // Računanje norme
    doublereal minimalna_norma;
    D = malloc(n_kvadrat * n * sizeof(doublereal));

    dgemm_("N", "N", &n_kvadrat, &n, &n, 
           &scalar_one, A, &n_kvadrat, Q, &n, 
           &scalar_zero, D, &n_kvadrat);

    for (integer i = 0; i < n_kvadrat * n; i++) {
        D[i] -= B[i];
    }

    minimalna_norma = dlange_("F", &n, &n_kvadrat, D, &n, WORK);
    printf("Minimalna norma je: %.5lf\n", minimalna_norma);


    return 0;
}
