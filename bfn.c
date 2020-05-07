#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bfn.h"

#define ZERO_THRESH 1e-12

bfn_t *bfn_new(int ang_mom, int nprim, double *e, double *c)
{
    bfn_t *bf;
    int n_nonzero_coef = 0;

    for (int i = 0; i < nprim; i++) {
        if (fabs(c[i]) > ZERO_THRESH) {
            n_nonzero_coef++;
        }
    }

    if (n_nonzero_coef == 0) {
        return NULL;
    }

    bf = (bfn_t *) malloc(sizeof(bfn_t));
    bf->e = (double *) malloc(sizeof(double) * n_nonzero_coef);
    bf->c = (double *) malloc(sizeof(double) * n_nonzero_coef);

    bf->l = ang_mom;
    bf->nprim = 0;
    for (int i = 0; i < nprim; i++) {
        if (fabs(c[i]) > ZERO_THRESH) {
            bf->e[bf->nprim] = e[i];
            bf->c[bf->nprim] = c[i];
            bf->nprim++;
        }
    }

    return bf;
}


void bfn_delete(bfn_t *bf)
{
    free(bf->e);
    free(bf->c);
    free(bf);
}


int bfn_same_exponents(bfn_t *b1, bfn_t *b2)
{
    if (b1->nprim != b2->nprim) {
        return 0;
    }

    for (int i = 0; i < b1->nprim; i++) {
        if (fabs(b1->e[i] - b2->e[i]) > ZERO_THRESH) {
            return 0;
        }
    }

    return 1;
}
