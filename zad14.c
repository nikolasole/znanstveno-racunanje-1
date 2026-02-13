#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

int main(int argc, char *argv[]) {
    // Parametri matrica i konstante
    integer broj_redaka = 4, broj_stupaca = 2, ione = 1, info;
    doublereal scalar_one = 1, scalar_zero = 0, scalar_mone = -1;

    // Alokacija memorije za matrice i vektore
    doublereal *matrica_A = malloc(broj_redaka * broj_stupaca * sizeof(doublereal));
    doublereal *matrica_B = malloc(broj_redaka * broj_stupaca * sizeof(doublereal));
    doublereal *matrica_Qa, *matrica_Qb, *matrica_C, *matrica_U, *matrica_Vt, *vektor_singularnih_vrijednosti;
    doublereal *vektor_tauA, *vektor_tauB, *radna_memorija, *radna_memorija_svd;
    doublereal *matrica_X, *matrica_Y;

    // Inicijalizacija matrica A i B
    matrica_A[0] = 1; matrica_A[1] = 0; matrica_A[2] = 0; matrica_A[3] = 0;
    matrica_A[4] = 1; matrica_A[5] = 1; matrica_A[6] = 1; matrica_A[7] = 1;

    matrica_B[0] = 1; matrica_B[1] = -1; matrica_B[2] = 1; matrica_B[3] = -1;
    matrica_B[4] = 0; matrica_B[5] = 1; matrica_B[6] = 0; matrica_B[7] = 1;

    // Alokacija memorije za QR faktorizaciju i SVD
    integer lwork = broj_stupaca * broj_stupaca;
    vektor_tauA = malloc(broj_stupaca * sizeof(doublereal));
    vektor_tauB = malloc(broj_stupaca * sizeof(doublereal));
    radna_memorija = malloc(lwork * sizeof(doublereal));

    matrica_C = malloc(broj_stupaca * broj_stupaca * sizeof(doublereal));
    vektor_singularnih_vrijednosti = malloc(broj_stupaca * sizeof(doublereal));
    matrica_U = malloc(broj_stupaca * broj_stupaca * sizeof(doublereal));
    matrica_Vt = malloc(broj_stupaca * broj_stupaca * sizeof(doublereal));

    // QR faktorizacija
    dgeqrf_(&broj_redaka, &broj_stupaca, matrica_A, &broj_redaka, vektor_tauA, radna_memorija, &lwork, &info);
    dgeqrf_(&broj_redaka, &broj_stupaca, matrica_B, &broj_redaka, vektor_tauB, radna_memorija, &lwork, &info);

    // Qa, Qb
    dorgqr_(&broj_redaka, &broj_stupaca, &broj_stupaca, matrica_A, &broj_redaka, vektor_tauA, radna_memorija, &lwork, &info);
    dorgqr_(&broj_redaka, &broj_stupaca, &broj_stupaca, matrica_B, &broj_redaka, vektor_tauB, radna_memorija, &lwork, &info);

    // C = Qa^T * Qb
    dgemm_("T", "N", &broj_stupaca, &broj_stupaca, &broj_redaka, 
           &scalar_one, matrica_A, &broj_redaka, matrica_B, &broj_redaka, 
           &scalar_zero, matrica_C, &broj_stupaca);

    // Singularna dekompozicija C
    lwork = 5 * broj_stupaca;
    radna_memorija_svd = malloc(lwork * sizeof(doublereal));

    dgesvd_("A", "A", &broj_stupaca, &broj_stupaca, matrica_C, &broj_stupaca, 
            vektor_singularnih_vrijednosti, matrica_U, &broj_stupaca, matrica_Vt, &broj_stupaca, 
            radna_memorija_svd, &lwork, &info);

    // Ispis
    printf("Kosinus glavnih kuteva:");
    for (integer i = 0; i < broj_stupaca; i++) {
        printf(" %.5lf", vektor_singularnih_vrijednosti[i]);
    }
    printf("\n");

    // X = A * U i Y = B * V
    matrica_X = malloc(broj_redaka * broj_stupaca * sizeof(doublereal));
    matrica_Y = malloc(broj_redaka * broj_stupaca * sizeof(doublereal));

    dgemm_("N", "N", &broj_redaka, &broj_stupaca, &broj_stupaca, 
           &scalar_one, matrica_A, &broj_redaka, matrica_U, &broj_stupaca, 
           &scalar_zero, matrica_X, &broj_redaka);

    dgemm_("N", "T", &broj_redaka, &broj_stupaca, &broj_stupaca, 
           &scalar_one, matrica_B, &broj_redaka, matrica_Vt, &broj_stupaca, 
           &scalar_zero, matrica_Y, &broj_redaka);

    // Ispis
    printf("\nMatrica X:\n");
    for (integer i = 0; i < broj_redaka; i++) {
        for (integer j = 0; j < broj_stupaca; j++) {
            printf("%.5lf ", matrica_X[i + broj_redaka * j]);
        }
        printf("\n");
    }

    printf("\nMatrica Y:\n");
    for (integer i = 0; i < broj_redaka; i++) {
        for (integer j = 0; j < broj_stupaca; j++) {
            printf("%.5lf ", matrica_Y[i + broj_redaka * j]);
        }
        printf("\n");
    }

    printf("Kut između ravnina iznosi: %.5lf radijana.\n", acos(vektor_singularnih_vrijednosti[1]));

    return 0;
}