#ifndef BASIS_LIB_H_INCLUDED
#define BASIS_LIB_H_INCLUDED

#include "basis.h"
#include "elements.h"

typedef struct {
    basis_t *basis_list[N_CHEM_ELEMENTS];
} basis_lib_t;

basis_lib_t *basis_lib_new();
void basis_lib_delete(basis_lib_t *);
basis_t *basis_lib_set(basis_lib_t *bas_lib, basis_t *bas);
basis_t *basis_lib_get(basis_lib_t *bas_lib, char *elem_symbol);
void basis_lib_print(basis_lib_t *bas_lib);

#endif /* BASIS_LIB_H_INCLUDED */
