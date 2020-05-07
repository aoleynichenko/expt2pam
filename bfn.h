#ifndef BFN_H_INCLUDED
#define BFN_H_INCLUDED

/**
 * contracted gaussian basis function
 */
typedef struct {
    int l;       /* angular momentum */
    int nprim;   /* length of expansion -- number of primitives */
    double *e;   /* exponents (alpha) */
    double *c;   /* contraction coefficients */
} bfn_t;

bfn_t *bfn_new(int ang_mom, int nprim, double *e, double *c);
void bfn_delete(bfn_t *bf);
int bfn_same_exponents(bfn_t *b1, bfn_t *b2);

#endif /* BFN_H_INCLUDED */
