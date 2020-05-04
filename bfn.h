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

#endif /* BFN_H_INCLUDED */
