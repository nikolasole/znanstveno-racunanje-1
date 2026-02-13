#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

main(integer argc, char *argv[])
{
    doublereal a[]={1.0,2.0,5.0,1.0, 1.0,1.0,1.0,4.0, 4.0,1.0,1.0,1.0, 1.0,6.0,0.0,3.0}; 
    doublereal b[]={7.0,10.0,7.0,9.0}; 
    doublereal x[]={1.0,1.0,1.0,1.0};
    doublereal alpha, rg;
    integer n, nhrs, *ipiv, info, incx;
    
    n=4;
    ipiv=malloc(n*sizeof(integer));
    nhrs=1;
    incx=1;
    alpha=-1.0;
    dgesv_( &n, &nhrs, a, &n, ipiv, b, &n, &info );
    printf("dgesv() je izracunao :[ %19.16f, %19.16f, %19.16f, %19.16f ]^T\n", b[0],b[1],b[2],b[3]);
    printf("Egzaktno rjesenje je : [ 1, 1, 1, 1 ]^T\n");
    daxpy_( &n,  &alpha, x, &incx, b, &incx );
    printf("Razlika izracunatog i egzaktnog rjesenja :[ %.2e, %.2e, %.2e, %.2e ]^T\n", b[0],b[1],b[2],b[3]);  
    rg=dnrm2_( &n, b, &incx )/dnrm2_( &n, x, &incx );
    printf("Relativna greska rjesenja iznosi : %.2e\n", rg);
}
