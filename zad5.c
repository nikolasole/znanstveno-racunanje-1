#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

// pomocne varijable
integer idist = 1;
integer iseed[] = {1,1,1,1};
integer INFO;
char sideL = 'L';
char sideR = 'R';
char transN = 'N';
char transT = 'T';
char TRANS = 'N';
doublereal scalar_one = 1;
integer one = 1;
doublereal zero = 0;
char uploL = 'L';
char uploU = 'U';
integer inc = 1;
char diagU = 'U';
char diagN = 'N';

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
    int k = 0;

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

// nekompletna faktorizacija Cholseskog
void ic (integer n, doublereal *a){
    int i,j,k;
    for (i=0; i<n; i++){
        for (k=0; k<i; k++)
            a[i+i*n]-=a[k+i*n]*a[k+i*n];
        a[i+i*n]=sqrt(a[i+i*n]);
        for (j=i+1; j<n; j++){
            if (a[i+j*n]!=0){
                for (k=0; k<i; k++)
                    a[i+j*n]-=a[k+i*n]*a[k+j*n];
                a[i+j*n]/=a[i+i*n];
            }        
        }
    }
}

// rjesava sustav Ax = b koristeći IC prekondicioniranu metodu konjugiranih gradijenata
void pcg(integer n, doublereal a[], doublereal b[], doublereal x[], doublereal tol){
    // alociranje memorije
    doublereal *m = (doublereal *)malloc(n*n*sizeof(doublereal));
    doublereal *r = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *d = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *p = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *ad = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *pom = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal rel_norma_rez;
    doublereal alpha;
    doublereal beta;
    doublereal stari=0;
    doublereal novi;
    doublereal norma_b;
    integer n_kvadrat = n*n;
    int k = 0;

    // kopiramo a u m
    dcopy_(&n_kvadrat, a, &inc, m, &inc);

    // racunamo normu od b
    norma_b = dnrm2_(&n, b, &inc);

    // r0 = b - A*x0 (x0 = 0, pa su jednaki b) // b->r
        dcopy_(&n, b, &inc, r, &inc);

    // riješi M*p0 = r0
        // vršimo IC
        ic(n, m);
        // kopiramo r u p
        dcopy_(&n, r, &inc, p, &inc);
        // rjesavamo s = R^(-T)*r
        dtrsv_(&uploU, &transT, &diagN, &n, m, &n, p, &inc);
        // rjesavamo p = R^(-1)*s
        dtrsv_(&uploU, &transN, &diagN, &n, m, &n, p, &inc);

    // d0 = p0 // p->d
        dcopy_(&n, p, &inc, d, &inc);

    // racunamo r*p (spremljen u stari)
    stari = ddot_(&n, r, &inc, p, &inc);

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

        // riješi M*p0 = r0
            // kopiramo r u p
            dcopy_(&n, r, &inc, p, &inc);
            // rjesavamo s = R^(-T)*r
            dtrsv_(&uploU, &transT, &diagN, &n, m, &n, p, &inc);
            // rjesavamo p = R^(-1)*s
            dtrsv_(&uploU, &transN, &diagN, &n, m, &n, p, &inc);

        // racunamo novi r*p
        novi = ddot_(&n, r, &inc, p, &inc);

        // racunamo beta = novi / stari
        beta = novi / stari;

        // racunamo d_k+1
            // kopiramo d u pom
            dcopy_(&n, d, &inc, pom, &inc);
            // kopiramo p u d
            dcopy_(&n, p, &inc, d, &inc);
            // d_k+1 = p_k+1(d) + beta*d(pom)
            daxpy_(&n, &beta, pom, &inc, d, &inc);

        // izjednačavamo stari = novi prije iduceg koraka
        stari = novi;

        // relativna norma reziduala (norma r / norma b)
        rel_norma_rez = dnrm2_(&n, r, &inc) / norma_b;

        k++;
        // u svakom koraku spisujemo relativnu normu reziduala
        printf("korak: %d, rel norma reziduala: %f \n", k, rel_norma_rez);
    }
    printf("x: ");
    for(int i=0;i<n;i++) printf("%f ",x[i]);
    printf("\n");
}

// rjesava sustav Ax = b koristeći dijagonalnu prekondicioniranu metodu konjugiranih gradijenata
void dpcg(integer n, doublereal a[], doublereal b[], doublereal x[], doublereal tol) {
    // alociranje memorije
    doublereal *r = (doublereal *)malloc(n * sizeof(doublereal));
    doublereal *d = (doublereal *)malloc(n * sizeof(doublereal));
    doublereal *p = (doublereal *)malloc(n * sizeof(doublereal));
    doublereal *ad = (doublereal *)malloc(n * sizeof(doublereal));
    doublereal *pom = (doublereal *)malloc(n * sizeof(doublereal));
    doublereal rel_norma_rez;
    doublereal alpha;
    doublereal beta;
    doublereal stari = 0;
    doublereal novi;
    doublereal norma_b;
    integer inc = 1;
    int k = 0;

    // racunamo normu od b
    norma_b = dnrm2_(&n, b, &inc);

    // r0 = b - A*x0 (x0 = 0, pa su jednaki b) // b->r
    dcopy_(&n, b, &inc, r, &inc);

    // Dijagonalno prekondicioniranje: p0 = r0 / A_ii
    for (int i = 0; i < n; i++) p[i] = r[i] / a[i + i * n];

    // d0 = p0 // p->d
    dcopy_(&n, p, &inc, d, &inc);

    // racunamo r*p (spremljen u stari)
    stari = ddot_(&n, r, &inc, p, &inc);

    // racunanje relativne norme reziduala (r=b, pa je 1)
    rel_norma_rez = 1;

    // algoritam
    while (rel_norma_rez > tol) {
        // racunamo i spemamo A*d u ad
        dsymv_(&sideL, &n, &scalar_one, a, &n, d, &inc, &zero, ad, &inc);

        // racunamo alpha_k = stari / d*A*d (d*ad)
        alpha = stari / ddot_(&n, d, &inc, ad, &inc);

        // racunamo x_k+1 = x_k + alpha*d
        daxpy_(&n, &alpha, d, &inc, x, &inc);

        // racunamo r_k+1
            // oduzimamo pa mijenjamo predznak od alpha
            alpha *= -1;
            // r_k+1 = r_k - alpha*ad
            daxpy_(&n, &alpha, ad, &inc, r, &inc);

        // Dijagonalno prekondicioniranje: p_k+1 = r_k+1 / A_ii
        for (int i = 0; i < n; i++) p[i] = r[i] / a[i + i * n];

        // racunamo novi r*p
        novi = ddot_(&n, r, &inc, p, &inc);

        // racunamo beta = novi / stari
        beta = novi / stari;

        // racunamo d_k+1
            // kopiramo d u pom
            dcopy_(&n, d, &inc, pom, &inc);
            // kopiramo p u d
            dcopy_(&n, p, &inc, d, &inc);
            // d_k+1 = p_k+1(d) + beta*d(pom)
            daxpy_(&n, &beta, pom, &inc, d, &inc);

        // izjednačavamo stari = novi prije iduceg koraka
        stari = novi;

        // relativna norma reziduala (norma r / norma b)
        rel_norma_rez = dnrm2_(&n, r, &inc) / norma_b;

        k++;
        // u svakom koraku spisujemo relativnu normu reziduala
        printf("korak: %d, rel norma reziduala: %f \n", k, rel_norma_rez);
    }
    printf("x: ");
    for (int i = 0; i < n; i++) printf("%f ", x[i]);
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

    // ucitavanje matrice A
    FILE *f;
    f=fopen("stieltjes_matr.txt","r");
    for(int i=0;i<n*n;i++) fscanf(f,"%lf",&a[i]);
    fclose(f);
        
    // racunamo b
        // x egzaktni = 1,1,1,...
        for(int i=0;i<n;i++) x_egzaktni[i]=1;

        // b = A*x_egzaktni
        dgemv_(&TRANS, &n, &n, &scalar_one, a, &n, x_egzaktni, &one, &zero, b, &one);

    // pozivamo metodu konjugiranih gradijenata
    //cg(n,a,b,x,tol);

    // pozivamo dijagonalno prekondicionirani sustav
    //dpcg(n,a,b,x,tol);

    // pozivamo IC prekondicioniranu metodu konjugiranih gradijenata
    pcg(n,a,b,x,tol);

    return 0;
}