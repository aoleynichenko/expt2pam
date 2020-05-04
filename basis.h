#ifndef BASIS_H_INCLUDED
#define BASIS_H_INCLUDED

#include "bfn.h"

typedef struct {
    int element;
    int nfun;
    bfn_t **functions;
} basis_t;

basis_t *basis_new(int element);
void basis_delete(basis_t *bas);
void basis_add_function(basis_t *bas, int ang_mom, int nprim, double *e, double *c);
void basis_print(basis_t *bas);

#endif /* BASIS_H_INCLUDED */
