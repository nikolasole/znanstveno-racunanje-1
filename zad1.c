#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

integer one = 1;
char uploL = 'L';
char uploR = 'U';
char side = 'L';
char transa = 'n';
char diag = 'n';
doublereal alpha = 1.0;
doublereal omega = 1.05;
doublereal epsilon = 1e-5;
int k = 0;
char norm = '1';
doublereal pom, normx0x, norma;
integer n_kvadrat;
doublereal mone = -1;

doublereal sor_norma(doublereal a[], doublereal b[], doublereal T[], doublereal M[], doublereal N[], doublereal c[], doublereal omega, integer n) {
    norma = 0.0;
    n_kvadrat = n * n;
    doublereal *WORK = (doublereal *)malloc(n * sizeof(doublereal));
    doublereal *L = (doublereal *)malloc(n * n * sizeof(doublereal));
    doublereal *R = (doublereal *)malloc(n * n * sizeof(doublereal));


    dlacpy_(&uploL, &n, &n, a, &n, L, &n);
    dlacpy_(&uploR, &n, &n, a, &n, R, &n);

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            if (i > j) {
                L[j * n + i] *= omega;
            }
        }
    }

    dcopy_(&n_kvadrat, L, &one, M, &one); 


    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            if (i < j) {
                R[j * n + i] *= omega * (-1);
            }
            if (i == j) {
                R[j * n + i] *= (1 - omega);
            }
        }
    }

    dcopy_(&n_kvadrat, R, &one, N, &one);


    dtrsm_(&side, &uploL, &transa, &diag, &n, &n, &alpha, L, &n, R, &n);
    dcopy_(&n_kvadrat, R, &one, T, &one);  


    dcopy_(&n, b, &one, c, &one);
    dtrsm_(&side, &uploL, &transa, &diag, &n, &n, &alpha, M, &n, c, &n);


    norma = dlange_(&norm, &n, &n, T, &n, WORK);

    return norma;
}

int sor_konvergencija(doublereal a[], doublereal b[], doublereal T[], doublereal M[], doublereal N[], doublereal c[], doublereal omega, integer n) {
    if (omega < 0 || omega >= 2) {
        return 0;
    } else if (sor_norma(a, b, T, M, N, c, omega, n) >= 1) {
        return 0;
    } else {
        return 1;
    }
}

void sor_rjesavac(doublereal a[], doublereal b[], doublereal x[], doublereal T[], doublereal M[], doublereal N[], doublereal c[], doublereal omega, doublereal epsilon, integer n, int *k) {
    doublereal *x0 = (doublereal *)malloc(n * sizeof(doublereal));
    doublereal *xpom = (doublereal *)malloc(n * sizeof(doublereal));
    doublereal *WORK = (doublereal *)malloc(n * sizeof(doublereal));

    dcopy_(&n, x, &one, x0, &one);

    for (int i = 0; i < n; i++) {
        x0[i] = (1 - omega) * x0[i];
        pom = b[i];

        for (int j = 0; j < i; j++)
            pom -= a[i + j * n] * x0[j];

        for (int j = i + 1; j < n; j++)
            pom -= a[i + j * n] * x0[j];

        x0[i] += pom * omega / a[i + i * n];
    }


    dcopy_(&n, x0, &one, xpom, &one);
    daxpy_(&n, &mone, x, &one, xpom, &one);
    normx0x = dlange_(&norm, &n, &one, xpom, &n, WORK);
    norma = sor_norma(a, b, T, M, N, c, omega, n);
    printf("norma T: %f \n", norma);
    *k = ceil(log((epsilon * (1 - norma)) / normx0x) / log(norma));

    for (int l = 1; l < *k; l++) {
        for (int i = 0; i < n; i++) {
            x0[i] = (1 - omega) * x0[i];
            pom = b[i];

            for (int j = 0; j < i; j++)
                pom -= a[i + j * n] * x0[j];
            for (int j = i + 1; j < n; j++)
                pom -= a[i + j * n] * x0[j];

            x0[i] += pom * omega / a[i + i * n];
        }
    }

    dcopy_(&n, x0, &one, x, &one);

}

int main() {
    integer n = 4;
    doublereal a[] = {101.0, -4.0, 8.0, 12.0, -4.0, 20.0, -7.0, 3.0, 8.0, -7.0, 78.0, 32.0, 12.0, 3.0, 32.0, 113.0};
    doublereal b[] = {117.0, 12.0, 111.0, 160.0};
    doublereal *T = (doublereal *)malloc(n * n * sizeof(doublereal));
    doublereal *M = (doublereal *)malloc(n * n * sizeof(doublereal));
    doublereal *N = (doublereal *)malloc(n * n * sizeof(doublereal));
    doublereal *c = (doublereal *)malloc(n * sizeof(doublereal));
    doublereal x_egzaktni[] = {1.0, 1.0, 1.0, 1.0};
    doublereal x[] = {0.0, 0.0, 0.0, 0.0};

    if (sor_konvergencija(a, b, T, M, N, c, omega, n)) {
        sor_rjesavac(a, b, x, T, M, N, c, omega, epsilon, n, &k);

        printf("Aproksimacija rješenja: [%f", x[0]);
        for (int i = 1; i < n; i++)
            printf(", %f", x[i]);
        printf("]\n");

        printf("Broj iteracija: %d \n", k);
    } else {
        printf("SOR ne konvergira \n");
    }

    return 0;
}