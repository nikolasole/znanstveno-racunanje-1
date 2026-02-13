#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

integer n=6;
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
doublereal omega = 1.35;
doublereal tol = 1e-8;
int k = 0;
char norm = 'I';
doublereal pom, normx0x, norma;
integer n_kvadrat;
doublereal mone = -1;
doublereal zero = 0;
doublereal kriterij;
doublereal A[]={11,-20,0,0,0,-2,-5,41,-3,0,-3,0,0,-15,7,-1,0,0,0,0,-4,2,-10,0,0,-6,0,-1,28,-15,-1,0,0,0,-15,47};
integer INFO;





static integer c__1 = 1;
static doublereal c_b7 = -1.;
static doublereal c_b8 = 1.;
static doublereal c_b20 = 0.;

int gmres_(n, b, x, restrt, work, ldw, h, ldh, iter, resid, matvec, psolve, 
           info)
   integer *n, *restrt, *ldw, *ldh, *iter, *info;
   doublereal *b, *x, *work, *h, *resid;
   int (*matvec) (), (*psolve) ();
{
    /* System generated locals */
    integer work_dim1, work_offset, h_dim1, h_offset, i__1;
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ int drot_();
    static doublereal bnrm2;
    extern doublereal dnrm2_();
    static integer i, k, r, s, v, w, y;
    extern /* Subroutine */ int dscal_(), basis_(), dcopy_(), drotg_();
    static integer maxit;
    static doublereal rnorm, aa, bb;
    static integer cs, av, sn;
    extern /* Subroutine */ int update_();
    static doublereal tol;

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = h_dim1 + 1;
    h -= h_offset;
    work_dim1 = *ldw;
    work_offset = work_dim1 + 1;
    work -= work_offset;
    --x;
    --b;

    /* Executable Statements */
    *info = 0;

/*     Test the input parameters. */

    if (*n < 0) {
	*info = -1;
    } else if (*ldw < max(1,*n)) {
	*info = -2;
    } else if (*iter <= 0) {
	*info = -3;
    } else if (*ldh < *restrt + 1) {
	*info = -4;
    }
    if (*info != 0) {
	return 0;
    }

    maxit = *iter;
    tol = *resid;

/*     Alias workspace columns. */

    r = 1;
    s = r + 1;
    w = s + 1;
    y = w;
    av = y;
    v = av + 1;

/*     Store the Givens parameters in matrix H. */

    cs = *restrt + 1;
    sn = cs + 1;

/*     Set initial residual (AV is temporary workspace here). */

    dcopy_(n, &b[1], &c__1, &work[av * work_dim1 + 1], &c__1);
    if (dnrm2_(n, &x[1], &c__1) != 0.) {

/*        AV is temporary workspace here. */

	dcopy_(n, &b[1], &c__1, &work[av * work_dim1 + 1], &c__1);
	(*matvec)(&c_b7, &x[1], &c_b8, &work[av * work_dim1 + 1]);
    }
    (*psolve)(&work[r * work_dim1 + 1], &work[av * work_dim1 + 1]);
    bnrm2 = dnrm2_(n, &b[1], &c__1);
    if (bnrm2 == 0.) {
	bnrm2 = 1.;
    }
    if (dnrm2_(n, &work[r * work_dim1 + 1], &c__1) / bnrm2 < tol) {
	goto L70;
    }

    *iter = 0;

L10:

    i = 0;

/*        Construct the first column of V. */

    dcopy_(n, &work[r * work_dim1 + 1], &c__1, &work[v * work_dim1 + 1], &
	    c__1);
    rnorm = dnrm2_(n, &work[v * work_dim1 + 1], &c__1);
    d__1 = 1. / rnorm;
    dscal_(n, &d__1, &work[v * work_dim1 + 1], &c__1);

/*        Initialize S to the elementary vector E1 scaled by RNORM. */

    work[s * work_dim1 + 1] = rnorm;
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	work[k + s * work_dim1] = 0.;
/* L20: */
    }

L30:

    ++i;
    ++(*iter);

    (*matvec)(&c_b8, &work[(v + i - 1) * work_dim1 + 1], &c_b20, &work[av * 
	    work_dim1 + 1]);
    (*psolve)(&work[w * work_dim1 + 1], &work[av * work_dim1 + 1]);

/*           Construct I-th column of H orthnormal to the previous */
/*           I-1 columns. */

    basis_(&i, n, &h[i * h_dim1 + 1], &work[v * work_dim1 + 1], ldw, &work[w *
	     work_dim1 + 1]);

/*           Apply Givens rotations to the I-th column of H. This */
/*           "updating" of the QR factorization effectively reduces */
/*           the Hessenberg matrix to upper triangular form during */
/*           the RESTRT iterations. */

    i__1 = i - 1;
    for (k = 1; k <= i__1; ++k) {
	drot_(&c__1, &h[k + i * h_dim1], ldh, &h[k + 1 + i * h_dim1], ldh, &h[
		k + cs * h_dim1], &h[k + sn * h_dim1]);
/* L40: */
    }

/*           Construct the I-th rotation matrix, and apply it to H so that
 */
/*           H(I+1,I) = 0. */

    aa = h[i + i * h_dim1];
    bb = h[i + 1 + i * h_dim1];
    drotg_(&aa, &bb, &h[i + cs * h_dim1], &h[i + sn * h_dim1]);
    drot_(&c__1, &h[i + i * h_dim1], ldh, &h[i + 1 + i * h_dim1], ldh, &h[i + 
	    cs * h_dim1], &h[i + sn * h_dim1]);

/*           Apply the I-th rotation matrix to [ S(I), S(I+1) ]'. This */
/*           gives an approximation of the residual norm. If less than */
/*           tolerance, update the approximation vector X and quit. */

    drot_(&c__1, &work[i + s * work_dim1], ldw, &work[i + 1 + s * work_dim1], 
	    ldw, &h[i + cs * h_dim1], &h[i + sn * h_dim1]);
    *resid = (d__1 = work[i + 1 + s * work_dim1], abs(d__1)) / bnrm2;
    printf("iteracija: %ld, rezidual: %lg\n", *iter, *resid); // ispis
    if (*resid <= tol) {
	update_(&i, n, &x[1], &h[h_offset], ldh, &work[y * work_dim1 + 1], &
		work[s * work_dim1 + 1], &work[v * work_dim1 + 1], ldw);
	goto L70;
    }
    if (*iter == maxit) {
	goto L50;
    }
    if (i < *restrt) {
	goto L30;
    }

L50:

/*        Compute current solution vector X. */

    update_(restrt, n, &x[1], &h[h_offset], ldh, &work[y * work_dim1 + 1], &
	    work[s * work_dim1 + 1], &work[v * work_dim1 + 1], ldw);

/*        Compute residual vector R, find norm, then check for tolerance. 
*/
/*        (AV is temporary workspace here.) */

    dcopy_(n, &b[1], &c__1, &work[av * work_dim1 + 1], &c__1);
    (*matvec)(&c_b7, &x[1], &c_b8, &work[av * work_dim1 + 1]);
    (*psolve)(&work[r * work_dim1 + 1], &work[av * work_dim1 + 1]);
    work[i + 1 + s * work_dim1] = dnrm2_(n, &work[r * work_dim1 + 1], &c__1);
    printf("iteracija: %ld, rezidual: %lg\n", *iter, *resid); // ispis
    *resid = work[i + 1 + s * work_dim1] / bnrm2;
    if (*resid <= tol) {
	goto L70;
    }
    if (*iter == maxit) {
	goto L60;
    }

/*        Restart. */

    goto L10;

L60:

/*     Iteration fails. */

    *info = 1;
    return 0;

L70:

/*     Iteration successful; return. */

    return 0;

/*     End of GMRES */

} /* gmres_ */


/*     =============================================================== */
/* Subroutine */ int update_(i, n, x, h, ldh, y, s, v, ldv)
integer *i, *n;
doublereal *x, *h;
integer *ldh;
doublereal *y, *s, *v;
integer *ldv;
{
    /* System generated locals */
    integer h_dim1, h_offset, v_dim1, v_offset;

    /* Local variables */
    extern /* Subroutine */ int dgemv_(), dcopy_(), dtrsv_();



/*     This routine updates the GMRES iterated solution approximation. */


/*     .. Executable Statements .. */

/*     Solve H*Y = S for upper triangualar H. */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    --s;
    --y;
    h_dim1 = *ldh;
    h_offset = h_dim1 + 1;
    h -= h_offset;
    --x;

    /* Function Body */
    dcopy_(i, &s[1], &c__1, &y[1], &c__1);
    dtrsv_("UPPER", "NOTRANS", "NONUNIT", i, &h[h_offset], ldh, &y[1], &c__1 );

/*     Compute current solution vector X = X + V*Y. */

    dgemv_("NOTRANS", n, i, &c_b8, &v[v_offset], ldv, &y[1], &c__1, &c_b8, &x[
	    1], &c__1 );

    return 0;

} /* update_ */


/*     ========================================================= */
/* Subroutine */ int basis_(i, n, h, v, ldv, w)
integer *i, *n;
doublereal *h, *v;
integer *ldv;
doublereal *w;
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1;
    doublereal d__1;

    /* Local variables */
    extern doublereal ddot_(), dnrm2_();
    static integer k;
    extern /* Subroutine */ int dscal_(), dcopy_(), daxpy_();



/*     Construct the I-th column of the upper Hessenberg matrix H */
/*     using the Gram-Schmidt process on V and W. */


    /* Parameter adjustments */
    --w;
    v_dim1 = *ldv;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    --h;

    /* Function Body */
    i__1 = *i;
    for (k = 1; k <= i__1; ++k) {
	h[k] = ddot_(n, &w[1], &c__1, &v[k * v_dim1 + 1], &c__1);
	d__1 = -h[k];
	daxpy_(n, &d__1, &v[k * v_dim1 + 1], &c__1, &w[1], &c__1);
/* L10: */
    }
    h[*i + 1] = dnrm2_(n, &w[1], &c__1);
    dcopy_(n, &w[1], &c__1, &v[(*i + 1) * v_dim1 + 1], &c__1);
    d__1 = 1. / h[*i + 1];
    dscal_(n, &d__1, &v[(*i + 1) * v_dim1 + 1], &c__1);

    return 0;
}


// y = α · A · x + β · y
int matvec(doublereal* alfa, doublereal* x, doublereal* beta, doublereal* y) //odi samo taj potprogram triba napraviti ali nez jel treba neke pointere tu smisljati ili ne
{
    dgemv_(&TRANS, &n, &n, alfa, A, &n, x, &one, beta ,y , &one);
}

// prekondicioniranje ne radimo
int psolve (doublereal* x, doublereal* b)
{
	dcopy_(&n, b, &one, x, &one);
    return 1;
}


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


int main() {
    // alociranje memorije
    doublereal *r = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *x = (doublereal *)malloc(n * sizeof(doublereal));
    integer restrt = n;
    integer ldh = restrt+1;
    integer iter = n;
    doublereal *h = (doublereal *)malloc((ldh*(restrt+2))*sizeof(doublereal));
    doublereal *WORK = (doublereal *)malloc((n)*(n+4)*sizeof(doublereal));

    // zadavanje vektora b
    doublereal b[] = {500, 0, 0, 0, 0, 0};

    // zadavanje x0 = 0
    for(int i=0;i<n;i++) x[i]=0;


    // ------ pozivamo metode ------
    //sor(A, b, x, r);
    //gs(A, b, x, r);
    gmres_(&n, b, x, &restrt, WORK, &n, h, &ldh, &iter, &tol, matvec, psolve, &INFO); 

    // ispis
    //printf("broj iteracija: %d \n", k);
    printf("Aproksimacija rješenja: [%f", x[0]);
        for (int i = 1; i < n; i++)
            printf(", %f", x[i]);
        printf("]\n");
    return 0;
}