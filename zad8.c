#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

integer n=99;
integer one = 1;
doublereal scalar_one = 1;
doublereal scalar_mone = -1;
char uploL = 'L';
char uploR = 'U';
char side = 'L';
char transa = 'n';
char TRANS = 'N';
char diag = 'n';
doublereal alpha = 1.0;
doublereal omega = 1.95;
doublereal tol = 1e-8;
int k = 0;
char norm = 'I';
doublereal pom, normx0x, norma;
integer n_kvadrat;
doublereal mone = -1;
doublereal zero = 0;
doublereal kriterij;
integer INFO;
integer inc = 1;


void sor(doublereal a[], doublereal b[], doublereal x[], doublereal r[]) {
    
    // racunamo normu od b koja se ponavlja
    doublereal norma_b = dnrm2_(&n, b, &one);

    // racunanje nultog kriterija
    dcopy_(&n, b, &one, r, &one);
    dgemv_(&TRANS, &n, &n, &scalar_mone, a, &n, x, &one, &scalar_one, r, &one);
    kriterij = dnrm2_(&n, r, &one)/norma_b;  

    while(kriterij>tol){
        for (int i=0;i<n;i++){
            x[i]=(1-omega)*x[i];
            pom=b[i];
            for (int j=0;j<i;j++)
                pom -= a[i+j*n]*x[j];
            for (int j=i+1;j<n;j++)
                pom-=a[i+j*n]*x[j];
            x[i]+=pom*omega/a[i+i*n];
        }
        k++;

        // racunanje kriterija
        dcopy_(&n, b, &one, r, &one);
        dgemv_(&TRANS, &n, &n, &scalar_mone, a, &n, x, &one, &scalar_one, r, &one);
        kriterij = dnrm2_(&n, r, &one)/norma_b;        
    }
}

void gs(doublereal a[], doublereal b[], doublereal x[], doublereal r[]) {
    
    // racunamo normu od b koja se ponavlja
    doublereal norma_b = dnrm2_(&n, b, &one);

    // racunanje nultog kriterija
    dcopy_(&n, b, &one, r, &one);
    dgemv_(&TRANS, &n, &n, &scalar_mone, a, &n, x, &one, &scalar_one, r, &one);
    kriterij = dnrm2_(&n, r, &one)/norma_b;  

    while(kriterij>tol){
        for (int i=0;i<n;i++){
            x[i]=b[i];
            for (int j=0;j<i;j++) x[i] = x[i]-a[i+j*n]*x[j];
            for (int j=i+1;j<n;j++) x[i] = x[i]-a[i+j*n]*x[j];
            x[i] = x[i]/a[i+i*n];
        }
        k++;

        // racunanje kriterija
        dcopy_(&n, b, &one, r, &one);
        dgemv_(&TRANS, &n, &n, &scalar_mone, a, &n, x, &one, &scalar_one, r, &one);
        kriterij = dnrm2_(&n, r, &one)/norma_b;        
    }
}

void cg(doublereal a[], doublereal b[], doublereal x[]){
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
    }
}


int main(int argc, char *argv[])
{
    // alociranje memorije
    doublereal *r = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *b = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *A = (doublereal *)malloc(n*n*sizeof(doublereal));
    doublereal *pom = (doublereal *)malloc(n*n*sizeof(doublereal));
    integer n_kvadrat = n*n;
    doublereal *y = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *y_egz = (doublereal *)malloc(n*sizeof(doublereal));

    // zadajemo y0
    for(int i=0; i<n; i++) y[i]=0;
    
    // zadajemo egzaktni y
    for(int i=0; i<n ; i++) y_egz[i]=0.01*(i+1)*cos(0.01*(i+1));

    // zadajemo A
    for(int i=0; i<n; i++){
	    if(i!=0) A[i*n+i-1]=-1;
	    A[i*n+i]=1.9999;
	    if(i!=n-1) A[i*n+i+1]=-1;
	}

    // zadajemo b
    for(int i=0; i<n; i++) b[i]=0.0002*sin(0.01*(i+1));
    b[98]+=cos(1);

    // provjeravamo pozitivnu definitnost
    for(int i=0; i<n_kvadrat; i++) pom[i]=A[i];
    dpotrf_(&uploL, &n, pom, &n, &INFO);
    if(INFO==0) printf("A je pozitivno definitna\n");
    else printf("A nije pozitivno definitna\n");


    // ------ pozivamo metode ------
    sor(A, b, y, r);
    //gs(A, b, y, r);
    //cg(A,b,y);
    
    // ispis
    printf("broj iteracija: %d \n", k);
    printf("Greška aproksimacije: [%f", y[0]-y_egz[0]);
        for (int i = 1; i < n; i++)
            printf(", %e", y[i]-y_egz[i]);
        printf("]\n");


    return 0;
}