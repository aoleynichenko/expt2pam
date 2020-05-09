#ifndef BASIS_H_INCLUDED
#define BASIS_H_INCLUDED

#include "elements.h"

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


typedef struct {
    int element;
    int nfun;
    bfn_t **functions;
} basis_t;

basis_t *basis_new(int element);
void basis_delete(basis_t *bas);
void basis_add_function(basis_t *bas, int ang_mom, int nprim, double *e, double *c);
void basis_print(basis_t *bas);


typedef struct {
    basis_t *basis_list[N_CHEM_ELEMENTS];
} basis_lib_t;

basis_lib_t *basis_lib_new();
void basis_lib_delete(basis_lib_t *);
basis_t *basis_lib_set(basis_lib_t *bas_lib, basis_t *bas);
basis_t *basis_lib_get(basis_lib_t *bas_lib, char *elem_symbol);
void basis_lib_print(basis_lib_t *bas_lib);

#endif /* BASIS_H_INCLUDED */
