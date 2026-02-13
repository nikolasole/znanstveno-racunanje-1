#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

int main(int argc, char *argv[]) {
    // Parametri matrica i konstante
    integer broj_redaka = 10, broj_stupaca = 2, ione = 1, info, radni_prostor_velicina;
    doublereal scalar_one = 1, scalar_zero = 0;

    // Alokacija memorije za matrice i vektore
    doublereal *matrica_A = malloc(broj_redaka * broj_stupaca * sizeof(doublereal));
    doublereal *vektor_Singularne = malloc(broj_stupaca * sizeof(doublereal));
    doublereal *matrica_U = malloc(broj_redaka * broj_redaka * sizeof(doublereal));
    doublereal *matrica_U_skracena = malloc(broj_redaka * broj_stupaca * sizeof(doublereal));
    doublereal *matrica_Vt = malloc(broj_stupaca * broj_stupaca * sizeof(doublereal));
    doublereal *vektor_radni_prostor;
    doublereal *vektor_pomocni = malloc(broj_stupaca * sizeof(doublereal));
    doublereal *vektor_koeficijenata = malloc(broj_stupaca * sizeof(doublereal));

    // Vrijednosti vektora b
    doublereal vektor_b[] = {3.5, 4.9, 6.8, 9.3, 10.9, 13.4, 15.1, 16.7, 19, 21.2};

    // Veličina radnog prostora za SVD
    radni_prostor_velicina = broj_redaka * 5;
    vektor_radni_prostor = malloc(radni_prostor_velicina * sizeof(doublereal));

    // Inicijalizacija matrice A
    for (integer i = 0; i < broj_redaka; i++) {
        matrica_A[i] = 1.0; // Prvi stupac matrice A
    }

    for (integer i = 0, j = 1; i < broj_redaka; i++, j++) {
        matrica_A[i + broj_redaka] = (doublereal)j; // Drugi stupac matrice A
    }

    // Singularna dekompozicija matrice A
    dgesvd_("S", "A", &broj_redaka, &broj_stupaca, matrica_A, &broj_redaka,
            vektor_Singularne, matrica_U, &broj_redaka, matrica_Vt, &broj_stupaca,
            vektor_radni_prostor, &radni_prostor_velicina, &info);

    // Ispis singularnih vrijednosti za provjeru
    printf("Singularne vrijednosti:\n");
    for (integer i = 0; i < broj_stupaca; i++) {
        printf("%.5lf\n", vektor_Singularne[i]);
    }

    // Računanje pom = U^T * b
    dgemv_("T", &broj_redaka, &broj_stupaca, &scalar_one, matrica_U, &broj_redaka,
           vektor_b, &ione, &scalar_zero, vektor_pomocni, &ione);

    // Normalizacija prema singularnim vrijednostima
    for (integer i = 0; i < broj_stupaca; i++) {
        vektor_pomocni[i] /= vektor_Singularne[i];
    }

    // Računanje koef = V * pom
    dgemv_("T", &broj_stupaca, &broj_stupaca, &scalar_one, matrica_Vt, &broj_stupaca,
           vektor_pomocni, &ione, &scalar_zero, vektor_koeficijenata, &ione);

    // Ispis aproksimativnog pravca
    printf("Aproksimativni pravac je: p(x) = %.5lf + %.5lf * x\n",
           vektor_koeficijenata[0], vektor_koeficijenata[1]);

    return 0;
}
