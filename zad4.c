#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>


// dlarnv
integer idist = 1;
integer iseed[] = {1,1,1,1};

// dgeqrf
integer INFO;

// dormqr
char sideL = 'L';
char sideR = 'R';
char transN = 'N';
char transT = 'T';

// dgemv
char TRANS = 'N';
doublereal scalar_one = 1;
integer one = 1;
doublereal zero = 0;

// dsymv
char uploL = 'L';

// ddot

// dscal

// daxpy

// dcopy
integer inc = 1;

// cg
int k = 0;


// rjesava sustav Ax = b koristeći metodu konjugiranih gradijenata
void cg(integer n, doublereal a[], doublereal b[], doublereal x[], doublereal tol){
    // alociranje memorije
    doublereal *r = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *d = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *ad = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *pom = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal rel_norma_rez;
    doublereal alpha;
    doublereal beta;
    doublereal stari=0;
    doublereal novi;
    doublereal norma_b;

    // racunamo normu od b
    norma_b = dnrm2_(&n, b, &inc);


    // d0 = r0 = b - A*x0 (x0 = 0, pa su jednaki b)
        // b->r
        dcopy_(&n, b, &inc, r, &inc);
        // b->d
        dcopy_(&n, b, &inc, d, &inc);

    // racunamo r*r (spremljen u stari)
    stari = ddot_(&n, r, &inc, r, &inc);

    // racunanje relativne norme reziduala (r=b, pa je 1)
    rel_norma_rez = 1;

    // algoritam
    while(rel_norma_rez > tol){
        // racunamo i spemamo A*d u ad
        dsymv_(&uploL, &n, &scalar_one, a, &n, d, &inc, &zero, ad, &inc);
        
        // racunamo alpha_k = stari / d*A*d (d*ad)
        alpha = stari / ddot_(&n, d, &inc, ad, &inc);

        // racunamo x_k+1 = x_k + alpha*d
        daxpy_(&n, &alpha, d, &inc, x, &inc);

        // racunamo r_k+1
            // oduzimamo pa mijenjamo predznak od alpha
            alpha*=-1;
            // r_k+1 = r_k - alpha*ad
            daxpy_(&n, &alpha, ad, &inc, r, &inc);

        // racunamo novi r*r
        novi = ddot_(&n, r, &inc, r, &inc);

        // racunamo beta = novi / stari
        beta = novi / stari;

        // racunamo d_k+1
            // kopiramo d u pom
            dcopy_(&n, d, &inc, pom, &inc);
            // kopiramo r u d
            dcopy_(&n, r, &inc, d, &inc);
            // d_k+1 = r_k+1(d) + beta*d(pom)
            daxpy_(&n, &beta, pom, &inc, d, &inc);

        // izjednačavamo stari = novi prije iduceg koraka
        stari = novi;

        // relativna norma reziduala (norma r = sqrt(r*r))
        rel_norma_rez = sqrt(novi) / norma_b;

        k++;
        // u svakom koraku spisujemo relativnu normu reziduala
        printf("korak: %d, rel norma reziduala: %f \n", k, rel_norma_rez);
    }
    printf("x: ");
    for(int i=0;i<n;i++) printf("%f ",x[i]);
    printf("\n");
}

int main(){
    // pocetni parametri
    integer n=100;
    integer n_kvadrat = n*n;
    doublereal tol = 1e-8;

    // alociranje memorije
    doublereal *a = (doublereal *)malloc(n*n*sizeof(doublereal));
    doublereal *b = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *x = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *x_egzaktni = (doublereal *)malloc(n*sizeof(doublereal));

    doublereal *q = (doublereal *)malloc(n*n*sizeof(doublereal));
    doublereal *WORK = (doublereal *)malloc(n*n*sizeof(doublereal));
    doublereal *TAU = (doublereal *)malloc(n*sizeof(doublereal));
    
    // pocetna iteracije x0 = 0,0,0,...
    for(int i=0;i<n;i++) x[i]=0;

    // zadavanje matrice A
        // zadavanje matrice LAMBDA (u A)
        for(int j=0;j<n;j++) for(int i=0;i<n;i++)
            if(i==j) a[j*n+i] = i/10+1;
            else a[j*n+i] = 0;
        // slucajna matrica (u Q)
        dlarnv_(&idist, iseed, &n_kvadrat, q);

        // slucajni ortogonalni Q (u Q)
        dgeqrf_(&n,&n,q,&n,TAU,WORK,&n_kvadrat,&INFO);

        // A = Q*LAMBDA*Qt
            // prvo slijeva
            dormqr_(&sideL, &transN, &n, &n, &n, q, &n, TAU, a, &n, WORK, &n_kvadrat, &INFO);
            // onda sdesna
            dormqr_(&sideR, &transT, &n, &n, &n, q, &n, TAU, a, &n, WORK, &n_kvadrat, &INFO);

    // racunamo b
        // x egzaktni = 1,1,1,...
        for(int i=0;i<n;i++) x_egzaktni[i]=1;

        // b = A*x_egzaktni
        dgemv_(&TRANS, &n, &n, &scalar_one, a, &n, x_egzaktni, &one, &zero, b, &one);

    // pozivamo metodu konjugiranih gradijenata
    cg(n,a,b,x,tol);

    return 0;
}