#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

integer n=100;
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
doublereal omega = 1.25;
doublereal epsilon = 1e-5;
int k = 0;
char norm = 'I';
doublereal pom, normx0x, norma;
integer n_kvadrat;
doublereal mone = -1;
doublereal zero = 0;
doublereal kriterij = 10;

void sor(doublereal a[], doublereal b[], doublereal x[], doublereal r[]) {
    
    // racunamo normu od b koja se ponavlja
    doublereal norma_b = dnrm2_(&n, b, &one);

    // racunanje nultog kriterija
    dcopy_(&n, b, &one, r, &one);
    dgemv_(&TRANS, &n, &n, &scalar_mone, a, &n, x, &one, &scalar_one, r, &one);
    kriterij = dnrm2_(&n, r, &one)/norma_b;  

    while(kriterij>epsilon){
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

int main() {
    // alociranje memorije
    doublereal *a = (doublereal *)malloc(n*n*sizeof(doublereal));
    doublereal *x_egzaktni = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *b = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *r = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *x = (doublereal *)malloc(n * sizeof(doublereal));
    doublereal *x_minus_x_egzaktni = (doublereal *)malloc(n * sizeof(doublereal));
    FILE *f;

    // ucitavanje matrice a iz fajla
    f=fopen("stieltjes_matr.txt","r");
    for(int i=0;i<n*n;i++) fscanf(f,"%lf",&a[i]);
    fclose(f);

    // zadavanje egzaktnog rjesenja x
    for(int i=0;i<n;i++) x_egzaktni[i]=1;

    // zadavanje vektora b
    dgemv_(&TRANS, &n, &n, &scalar_one, a, &n, x_egzaktni, &one, &zero, b, &one);

    // zadavanje x0 = 0
    for(int i=0;i<n;i++) x[i]=0;

    sor(a, b, x, r);

    // racunanje x - x_egzaktni
    for(int i=0;i<n;i++) x_minus_x_egzaktni[i]=x[i]-x_egzaktni[i];

    printf("relativna greška: %e \n", dnrm2_(&n, x_minus_x_egzaktni, &one)/dnrm2_(&n, x_egzaktni, &one));
    printf("Broj iteracija: %d \n", k);

    return 0;
}