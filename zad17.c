#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

doublereal f(doublereal x, doublereal y){
    return (x*x*y - x*x - y*y + 175)/250;
}

int main(void){

    integer n = 21, rang, ione=1;
    doublereal x, y, nul, jed, minjed;
    int k;
    doublereal *A = malloc(n*n*sizeof(doublereal));
    doublereal *Q = malloc(n*n*sizeof(doublereal));
    doublereal *b = malloc(n*sizeof(doublereal));
    doublereal *x0= malloc(n*sizeof(doublereal));
    doublereal *x1= malloc(n*sizeof(doublereal));
    integer *jpvt = malloc(n*sizeof(int));
    doublereal *tau = malloc(n*sizeof(doublereal));
    integer lwork = 3*n+1, info = 0;
    doublereal *work = malloc(lwork*sizeof(doublereal));
    doublereal tol = 21*1e-16;
    doublereal *diag = malloc(n*sizeof(doublereal));
    doublereal min;

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            x = -5 + 0.5*i;
            y = -5 + 0.5*j;
            A[i + j*n] = f(x,y);
        }
    }

    printf("Matrica A: \n");
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            printf("%lf    ", A[i+j*n]);
        }
        printf("\n");
    }
    printf("\n\n");

    dlaset_("a", &n, &ione, &jed, &jed, b, &n);
    b[0] = 10;
    dlacpy_("a", &n, &n, A, &n, Q, &n);
    dgeqp3_(&n, &n, Q, &n, jpvt, tau, work, &lwork, &info);
    dlacpy_("a", &n, &n, Q, &n, A, &n);
    dorgqr_(&n, &n, &n, Q, &n, tau, work, &lwork, &info);

    printf("Matrica jpvt: \n");
    for(int i=0; i<n; i++){
        printf("%ld  ", jpvt[i]);
    }
    printf("\n\n");

    for (int i=0; i<n; i++){
        if (abs(A[i+n*i]) <= tol) {
            rang = i;
            break;
        }
    }
    printf("rang = %ld\n\n", rang);

    doublereal *R = malloc(rang*n*sizeof(doublereal));
    doublereal *Z = malloc(rang*n*sizeof(doublereal));

    for(int i=0;i<rang; i++){
        for(int j=i; j<n; j++){
            R[i+j*rang] = A[i+j*n];
        }
    }

    printf("Matrica A nakon QR: \n");
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            printf("%lg  ", A[i+j*n]);
        }
        printf("\n");
    }
    printf("\n\n");

    printf("Matrica R: \n");
    for(int i=0; i<rang; i++){
        for(int j=0; j<n; j++){
            printf("%lf ", R[i+j*rang]);
        }
        printf("\n");
    }
    printf("\n\n");

    printf("Matrica Q: \n");
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            printf("%lg  ", Q[i+j*n]);
        }
        printf("\n");
    }
    printf("\n\n");

    doublereal *tau1 = malloc(rang*sizeof(doublereal));

    dlacpy_("a", &rang, &n, R, &rang, Z, &rang);
    dgelqf_(&rang, &n, Z, &rang, tau1, work, &lwork, &info);
    dlacpy_("a", &rang, &n, Z, &rang, R, &rang);
    dorglq_(&rang, &n, &rang, Z, &rang, tau1, work, &lwork, &info);

    doublereal *L = malloc(rang*rang*sizeof(doublereal));
    doublereal *QR = malloc(n*rang*sizeof(doublereal));
    doublereal *pom = malloc(rang*sizeof(doublereal));
    doublereal *pom2 = malloc(rang*sizeof(doublereal));

    printf("Matrica Z: \n");
    for(int i=0; i<rang; i++){
        for(int j=0; j<n; j++){
            printf("%lg  ", Z[i+j*rang]);
        }
        printf("\n");
    }
    printf("\n\n");

    printf("Matrica R nakon LQ: \n");
    for(int i=0; i<rang; i++){
        for(int j=0; j<n; j++){
            printf("%lg  ", R[i+j*rang]);
        }
        printf("\n");
    }
    printf("\n\n");

    for(int i=0;i<rang; i++){
        for(int j=0; j<i+1; j++){
            L[i+j*rang] = R[i+j*rang];
        }
    }

    printf("Matrica L11:\n");
    for(int i=0; i<rang; i++) {
        for(int j=0;j<rang; j++){
            printf("%lf  ", L[i+j*rang]);
        }
        printf("\n");
    }
    printf("\n\n");

    for(int i=0; i<n; i++){
        for(int j=0; j<rang; j++){
            QR[i+j*n] = Q[i + j*n];
        }
    }

    printf("Matrica Q(:,1:r):\n");
    for(int i=0; i<n; i++) {
        for(int j=0;j<rang; j++){
            printf("%lf  ", Q[i+j*n]);
        }
        printf("\n");
    }
    printf("\n\n");

    dgemv_("T", &n, &rang, &jed, QR, &n, b, &ione, &nul, pom, &ione);

    printf("pom nakon prvog množenja:\n");
    for(int i=0; i<rang; i++){
        printf("%lf\n", pom[i]);
    }
    printf("\n\n");

    dtrsm_("l", "l", "n", "n", &rang, &ione, &jed, L, &rang, pom, &rang);

    printf("pom nakon drugog množenja:\n");
    for(int i=0; i<rang; i++){
        printf("%lf\n", pom[i]);
    }
    printf("\n\n");

    dgemv_("T", &rang, &n, &jed, Z, &rang, pom, &ione, &nul, x0, &ione);

    printf("Rješenje x prije pivotiranja:\n");
    for(int i=0; i<n; i++){
        printf("%lf\n", x0[i]);
    }
    printf("\n\n");

    dlacpy_("a", &n, &ione, x0, &n, x1, &n);

    for(int i=0; i<n; i++){
        k = jpvt[i];
        x1[k] = x0[i];
    }

    printf("Rješenje x:\n");
    for(int i=0; i<n; i++){
        printf("%lf\n", x1[i]);
    }
    printf("\n\n");
}
