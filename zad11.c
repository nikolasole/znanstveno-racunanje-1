#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

int main(int argc, char *argv[]) {
    // Dimenzije matrica i parametri
    integer n = 4;
    integer n_kvadrat = 16;
    integer INFO;
    integer inc = 1;

    doublereal scalar_mone = -1;
    doublereal scalar_one = 1;
    doublereal zero = 0;

    // Alokacija memorije
    doublereal *A = malloc(n_kvadrat * sizeof(doublereal));
    doublereal *Apom = malloc(n_kvadrat * sizeof(doublereal));
    doublereal *x = malloc(n_kvadrat * sizeof(doublereal));
    doublereal *w = malloc(n * sizeof(doublereal));
    doublereal *WORK = malloc(n_kvadrat * sizeof(doublereal));
    doublereal *Mminp = malloc(n_kvadrat * sizeof(doublereal));

    // Inicijalizacija matrica M i K
    doublereal M[] = {2, 0, 0, 0, 
                      0, 5, 0, 0, 
                      0, 0, 3, 0, 
                      0, 0, 0, 6};

    doublereal K[] = {24, -9, -5,  0,
                      -9, 22, -8, -5,
                      -5, -8, 25, -7,
                       0, -5, -7, 18};

    // Izračunavanje M^(-1/2)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if ((i + n * j) % (n + 1) == 0) {
                Mminp[i + n * j] = 1 / sqrt(M[i + n * j]);
            } else {
                Mminp[i + n * j] = 0;
            }
        }
    }

    // Izračunavanje A = M^(-1/2) * K * M^(-1/2)
    dgemm_("N", "N", &n, &n, &n, &scalar_one, Mminp, &n, K, &n, &zero, Apom, &n);
    dgemm_("N", "N", &n, &n, &n, &scalar_one, Apom, &n, Mminp, &n, &zero, A, &n);

    // Provjera matrice A
    printf("Matrica A za provjeru:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf(" %8.4lf ", A[i + j * n]);
        }
        printf("\n");
    }

    // Računanje svojstvenih vrijednosti i svojstvenih vektora
    dsyev_("V", "U", &n, A, &n, w, WORK, &n_kvadrat, &INFO);

    printf("\nSvojstvene vrijednosti i vektori:\n");
    for (int i = 0; i < n; i++) {
        printf("Svojstvena vrijednost %d: %.4lf\n", i + 1, w[i]);
    }
    printf("\n");

    // Računanje konačnog rješenja x
    dgemm_("N", "N", &n, &n, &n, &scalar_one, Mminp, &n, A, &n, &zero, x, &n);
    for (int i = 0; i < n; i++) {
        printf("x%d = [", i + 1);
        for (int j = 0; j < n; j++) {
            printf(" %.4lf", x[j + i * n]);
        }
        printf("] * e^(%.4lf * it)\n", sqrt(w[i]));
    }

    return 0;
}
