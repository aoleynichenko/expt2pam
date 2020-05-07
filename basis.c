#include <stdio.h>
#include <stdlib.h>

#include "basis.h"
#include "bfn.h"

basis_t *basis_new(int element)
{
    basis_t *bas;

    bas = (basis_t *) malloc(sizeof(basis_t));

    bas->element = element;
    bas->nfun = 0;
    bas->functions = NULL;

    return bas;
}


void basis_delete(basis_t *bas)
{
    for (int i = 0; i < bas->nfun; i++) {
        bfn_t *bf = bas->functions[i];
        bfn_delete(bf);
    }
    free(bas->functions);
    free(bas);
}


void basis_add_function(basis_t *bas, int ang_mom, int nprim, double *e, double *c)
{
    bfn_t *bf = bfn_new(ang_mom, nprim, e, c);
    if (bf == NULL) {
        return;
    }

    if (bas->nfun == 0) {
        bas->nfun = 1;
        bas->functions = (bfn_t **) malloc(sizeof(bfn_t *) * 1);
    }
    else {
        bas->nfun++;
        bas->functions = (bfn_t **) realloc(bas->functions, sizeof(bfn_t *) * bas->nfun);
    }

    bas->functions[bas->nfun - 1] = bf;
}


void basis_print(basis_t *bas)
{
    printf("Basis set for E%d\n", bas->element);
    printf("Number of contracted basis functions: %d\n", bas->nfun);
    for (int i = 0; i < bas->nfun; i++) {
        bfn_t *bf = bas->functions[i];
        printf("[%d] L=%d nprim=%d\n", i, bf->l, bf->nprim);
        for (int j = 0; j < bf->nprim; j++) {
            printf("\t%18.8e%18.8e\n", bf->e[j], bf->c[j]);
        }
    }
}
