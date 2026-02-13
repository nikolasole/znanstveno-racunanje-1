#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

integer n=100;
doublereal *A;
integer INFO;

// dlarnv
integer idist = 1;
integer iseed[] = {1,1,1,1};

// dgemv
char TRANS = 'N';
doublereal scalar_one = 1;
integer one = 1;
doublereal zero = 0;

// dgeev
char JOBVL = 'N';
char JOBVR = 'N';


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

int main(){
    // pocetni parametri
    integer n_kvadrat = n*n;
    doublereal tol = 1e-5;

    // alociranje memorije
    A=malloc(n_kvadrat*sizeof(doublereal));
    doublereal *x = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *f = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *g = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *b = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *V = (doublereal *)malloc(n_kvadrat*sizeof(doublereal));
    doublereal *B = (doublereal *)malloc(n_kvadrat*sizeof(doublereal));
    doublereal *pom = (doublereal *)malloc(n_kvadrat*sizeof(doublereal));
    doublereal *Bpom = (doublereal *)malloc(n_kvadrat*sizeof(doublereal));
    doublereal *TAU = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *WORK = (doublereal *)malloc(n*n*sizeof(doublereal));
    integer *IPIV = (integer *)malloc(n*sizeof(integer));
    
    integer restrt = n;
    integer ldh = restrt+1;
    integer iter = n;
    doublereal *h = (doublereal *)malloc((ldh*(restrt+2))*sizeof(doublereal));

    doublereal *VL = (doublereal *)malloc(n_kvadrat*sizeof(doublereal));
    doublereal *VR = (doublereal *)malloc(n_kvadrat*sizeof(doublereal));
    doublereal *WR = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *WI = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *arg = (doublereal *)malloc(n*sizeof(doublereal));
    doublereal *z = (doublereal *)malloc(n*sizeof(doublereal));


    // postavljamo funkciju f
    for(int i=0; i<n; i++) f[i]=(n-i);

    // postavljamo funkciju g
    for(int i=1; i<n+1 ;i++) g[i-1]=sqrt((n-i+1)*(n-i+1)-(n-i)*(n-i));

    // zadavanje V
        // slucajna matrica (u V)
        dlarnv_(&idist, iseed, &n_kvadrat, V);

        // slucajni ortogonalni V (u V)
        dgeqrf_(&n, &n, V, &n, TAU, WORK, &n_kvadrat, &INFO);

        // spaja ovaj TAU u V?
        dorgqr_(&n, &n, &n, V, &n, TAU, WORK, &n_kvadrat, &INFO);

    // b=V*g 
    dgemv_(&TRANS, &n, &n, &scalar_one, V, &n, g, &one, &zero, b, &one);
    
    // B = [b v1 · · · v99]
    for (int i=0; i<n_kvadrat; i++) if (i<n) B[i]=b[i]; else B[i]=V[i-n];

    // zadavanje A
        // pomocna matrica u sredini; pom; jedinice na poddijagonali i gore desno
        for (int j=0; j<n; j++) for (int i=0; i<n; i++){
            pom[i+j*n]=0; 
            if (i==j+1) pom[n*j+i]=1;
            if (i==0 && j==n-1) pom[n*j+i]=1;
        }
    
        // A=B*pom*B^-1
            // Bpom = B*pom
            dgemm_(&TRANS, &TRANS, &n, &n, &n, &scalar_one, B, &n, pom, &n, &zero, Bpom, &n);
        
            // racunamo B^-1, spremljen u B
            dgetrf_(&n ,&n, B, &n, IPIV, &INFO); 
            dgetri_(&n, B, &n, IPIV, WORK, &n_kvadrat, &INFO);
            
            // A = Bpom*B^-1
            dgemm_(&TRANS, &TRANS, &n, &n, &n, &scalar_one, Bpom, &n, B, &n, &zero, A, &n);
        
    // postavljamo x na 0
    for(int i=0; i<n; i++) x[i]=0;

    // pozivamo GMRES
    gmres_(&n, b, x, &restrt, WORK, &n, h, &ldh, &iter, &tol, matvec, psolve, &INFO);

    // provjeravamo svojstvene vrijednosti
        // racunamo svojstvene vrijednosti od A, kompleksne, spremljene u WR i WI
        dgeev_(&JOBVL,&JOBVR, &n, A, &n, WR, WI, VL, &n, VR, &n, WORK, &n_kvadrat, &INFO);
    
        for(int i=0; i<n; i++){
            // prebacujemo u polarni oblik
            arg[i] = atan2(WI[i], WR[i]);
            z[i] = sqrt(WR[i]*WR[i] + (WI[i]*WI[i]));

            // dizemo na 100tu potenciju
            arg[i] = arg[i]*100; 
            z[i] = pow(z[i],100);

            // usporedujemo s jedinicama?
            printf("%lg + %lg*i \n", z[i]*cos(arg[i]), z[i]*sin(arg[i]));
        }

    return 0;
}