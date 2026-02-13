#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

integer n;
integer iseed[] = {1,1,1,1};
integer one = 1;
doublereal scalar_one = 1;
doublereal scalar_mone = -1;
char uploL = 'L';
char uploR = 'U';
char uploN = 'N';
char sideL = 'L';
char sideR = 'R';
char transN = 'N';
char transT = 'T';
char diag = 'n';
doublereal alpha;
doublereal beta;
doublereal tol;
int k = 0;
char norm = 'I';
doublereal pom, normx0x, norma;
integer n_kvadrat;
doublereal mone = -1;
doublereal zero = 0;
doublereal kriterij;
integer INFO;
integer inc = 1;
integer idist=3;


doublereal sgn(doublereal x){
    if (x > 0) return 1;
    else if (x < 0) return -1;
    else return 0;
}

void jacobi_sd(integer n, doublereal *A, doublereal tol){
    doublereal app, apq, aqq, TAU, t, *work, s, c, pom;
    doublereal normA, S_A, Spom=0;
    work=malloc(n*sizeof(doublereal));
    normA=dlange_( "F", &n, &n, A, &n, work);
    int iter=0;
    
    //racunanje S(A) vandijgonalna norma
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            if((i+j*n)%(n+1)!=0){
                Spom+=A[i+j*n]*A[i+j*n];
            }
        } 
    }
    S_A=sqrt(Spom);
    printf("pomocna suma vandijagonalnih elemenata: %lf \n", S_A);
    while (S_A > tol*normA ){
        iter++;
        for (int p = 0; p < n - 1; p++){
            for (int q = p + 1; q < n; q++){
                if (A[p+q*n] != 0){
                    TAU = (A[q+q*n] - A[p+p*n])/(2 * A[p+q*n]);
                    t =sgn(TAU)/(abs(TAU) + sqrt(1 + pow(TAU, 2)));
                    c = 1/(sqrt(1 + pow(t, 2)));
                    s = t * c;
                    app = A[p+n*p]; apq = A[p+n*q]; aqq = A[q+n*q];
                    app = app - t * apq;
                    aqq = aqq + t * apq;
                    for (int k = 0; k < n; k++){
                        pom = A[k+n*p];
                        A[k+n*p] = c * pom - s * A[k+n*q];
                        A[k+n*q] = s * pom + c * A[k+n*q];
                        A[p+n*k] = A[k+n*p]; A[q+n*k] = A[k+n*q];
                        }
                    A[p+n*q] = 0; A[q+n*p] = 0;
                    A[p+n*p] = app; A[q+n*q] = aqq;
                }
            }
        }
	
        //racunanje S(A) vandijgonalna norma
        Spom=0;
        for(int i=0; i<n; i++)
        {
            for(int j=0; j<n; j++)
            {
                if((i+j*n)%(n+1)!=0){
                    Spom+=A[i+j*n]*A[i+j*n];
                }
            } 
        }
    S_A=sqrt(Spom);
    printf("pomocna suma vandijagonalnih elemenata: %.15lf \n", S_A);
        
    }
    printf("broj iteracija: %d \n", iter);
}

int main(int argc, char *argv[])
{
    // a
        printf("Podzadatak a: \n");

        // alociranje memorije
        n=4; n_kvadrat=n*n;
        tol = 4*1e-16;
        doublereal *r = (doublereal *)malloc(n*sizeof(doublereal));
        doublereal *b = (doublereal *)malloc(n*sizeof(doublereal));
        doublereal *A = (doublereal *)malloc(n*n*sizeof(doublereal));
        doublereal *SlucU = (doublereal *)malloc(n*n*sizeof(doublereal));
        doublereal *x = (doublereal *)malloc(n*sizeof(doublereal));
        doublereal *x0 = (doublereal *)malloc(n*sizeof(doublereal));
        doublereal *TAU = (doublereal *)malloc(n*sizeof(doublereal));
        doublereal *WORK = (doublereal *)malloc(n*sizeof(doublereal));

        // zadajemo A
        alpha = 0;
        beta = 1;
        dlaset_(&uploN, &n, &n, &alpha, &beta, A, &n);
        A[0]=-10;
        A[5]=-5;
        A[10]=0.1;
        A[15]=0.2;

        dlarnv_(&idist, iseed, &n_kvadrat, SlucU);

        dgeqrf_(&n, &n, SlucU, &n, TAU, WORK, &n, &INFO); 
        dormqr_(&sideL, &transN, &n, &n, &n, SlucU, &n, TAU, A, &n, WORK, &n, &INFO);
        dormqr_(&sideR, &transT, &n, &n, &n, SlucU, &n, TAU, A, &n, WORK, &n, &INFO);
	
        for(int i=0; i<n; i++){
            x[i]=1;
            x0[i]=0;
        }

        alpha = 1;
        beta = 0;
        dgemv_(&transN, &n, &n, &alpha, A, &n, x, &one, &beta, b, &one);

	// pokrecemo jacobi sd
        jacobi_sd(n, A, tol);

        // ispis
        printf("Dijagonalni elementi: %lf %lf %lf %lf \n",A[0], A[5], A[10], A[15]);


	
    // b
	printf("\nPodzadatak b: \n");

	// alociranje memorije
	integer n1=10;
	doublereal tol1=1e-15;
	doublereal *Ries = (doublereal *)malloc(n1*n1*sizeof(doublereal));

	// zadajemo Riesovu matricu
	for(int i=0; i<n1; i++) for(int j=0; j<n1; j++) Ries[i+j*n1]=1/(2*(8-i-j+1.5));
	
	// pokrecemo jacobi sd 
	jacobi_sd(n1, Ries, tol1);

	// ispis
	printf("Dijagonalni elementi: ");
	for(int i=0; i<n1; i++) printf("%lf  ", Ries[i+i*n1]);
    

    return 0;
}
