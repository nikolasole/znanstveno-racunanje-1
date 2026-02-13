#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

int main(int argc, char *argv[]) {
    // Definicije parametara
    integer n = 7;
    integer ione = 1;
    integer dva = 2;
    integer *v1, *v2, *iwork, *ifail;
    
    doublereal nul = 0;
    doublereal scalar_one = 1;
    doublereal scalar_mone = -1;
    doublereal *W, *U, *Z, *D, *L, *LN, *A, *w, *WORK, abstol;

    // Alokacija memorije
    W = malloc(n * n * sizeof(doublereal));
    U = malloc(n * sizeof(doublereal));
    v1 = malloc(n * sizeof(integer));
    v2 = malloc(n * sizeof(integer));
    Z = malloc(n * n * sizeof(doublereal));
    D = malloc(n * n * sizeof(doublereal));
    L = malloc(n * n * sizeof(doublereal));
    LN = malloc(n * n * sizeof(doublereal));
    A = malloc(n * n * sizeof(doublereal));
    w = malloc(n * sizeof(doublereal));
    WORK = malloc(8 * n * sizeof(doublereal));
    iwork = malloc(5 * n * sizeof(integer));
    ifail = malloc(n * sizeof(integer));

    // Inicijalizacija tolerancije
    abstol = 2 * dlamch_("S");

    integer M, info, n_kvadrat = 8 * n;

    // Postavljanje matrice W
    dlaset_("A", &n, &n, &nul, &nul, W, &n);
    W[1] = 2;
    W[2] = 3;
    W[3] = 4;
    W[7] = 2;
    W[10] = 7;
    W[11] = 1;
    W[14] = 3;
    W[17] = 3;
    W[19] = 2;
    W[20] = 1;
    W[21] = 4;
    W[22] = 7;
    W[23] = 3;
    W[29] = 1;
    W[33] = 7;
    W[34] = 3;
    W[37] = 2;
    W[39] = 7;
    W[41] = 5;
    W[44] = 1;
    W[46] = 3;
    W[47] = 5;

    // Postavljanje vektora w koji čuva težine vrhova
    for (int i = 0; i < n; i++) {
        w[i] = 0;
        for (int j = 0; j < n; j++) {
            w[i] += W[i + j * n];
        }
    }

    // Komadi za diag
    dlaset_("A", &n, &n, &nul, &nul, D, &n);
    for (int i = 0; i < n; i++) {
        D[i + n * i] = w[i];
    }

    // Laplaceova matrica grafa
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            L[i + j * n] = D[i + j * n] - W[i + j * n];
        }
    }

    // Računanje D^-1/2
    for (int i = 0; i < n; i++) {
        D[i + n * i] = 1 / sqrt(D[i + n * i]);
    }

    // Računanje normalizirane LN
    dgemm_("N", "N", &n, &n, &n, &scalar_one, D, &n, L, &n, &nul, A, &n);
    dgemm_("N", "N", &n, &n, &n, &scalar_one, A, &n, D, &n, &nul, LN, &n);

    // Prvi dsyevx poziv
    dsyevx_("V", "I", "U", &n, L, &n, &scalar_one, &scalar_one, &dva, &dva, &abstol, &M, U, Z, &n, WORK, &n_kvadrat, iwork, ifail, &info);
    printf("Svojstvena vrijednost: L2 = %.6lf\n", U[0]);

    printf("Fiedlerov vektor u2:\n");
    for (int i = 0; i < n; i++) {
        printf("%.6lf\n", Z[i]);
    }

    // Klasifikacija vrhova prema Fiedlerovom vektoru
    int j = 0, k = 0;
    for (int i = 0; i < n; i++) {
        if (Z[i] >= 0) {
            v1[j++] = i + 1;
        } else {
            v2[k++] = i + 1;
        }
    }

    printf("V1 = {");
    for (int i = 0; i < j; i++) {
        printf("%ld%s", v1[i], (i < j - 1) ? ", " : "}\n");
    }

    printf("V2 = {");
    for (int i = 0; i < k; i++) {
        printf("%ld%s", v2[i], (i < k - 1) ? ", " : "}\n");
    }

    // Drugi dsyevx poziv
    dsyevx_("V", "I", "U", &n, LN, &n, &scalar_one, &scalar_one, &dva, &dva, &abstol, &M, U, Z, &n, WORK, &n_kvadrat, iwork, ifail, &info);
    printf("\nSvojstvena vrijednost: LN2 = %.6lf\n", U[0]);

    printf("Fiedlerov vektor uN2:\n");
    for (int i = 0; i < n; i++) {
        printf("%.6lf\n", Z[i]);
    }

    for (int i = 0; i < n; i++) {
        Z[i] += Z[i] * D[i + i * n];
    }

    j = 0;
    k = 0;
    for (int i = 0; i < n; i++) {
        if (Z[i] >= 0) {
            v1[j++] = i + 1;
        } else {
            v2[k++] = i + 1;
        }
    }

    printf("V1 = {");
    for (int i = 0; i < j; i++) {
        printf("%ld%s", v1[i], (i < j - 1) ? ", " : "}\n");
    }

    printf("V2 = {");
    for (int i = 0; i < k; i++) {
        printf("%ld%s", v2[i], (i < k - 1) ? ", " : "}\n");
    }

    return 0;
}
