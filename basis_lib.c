#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basis.h"
#include "basis_lib.h"
#include "error.h"

basis_lib_t *basis_lib_new()
{
    basis_lib_t *bas_lib = (basis_lib_t *) malloc(sizeof(basis_lib_t));

    memset(bas_lib->basis_list, 0, sizeof(bas_lib->basis_list));

    return bas_lib;
}


void basis_lib_delete(basis_lib_t *bas_lib)
{
    int sz = sizeof(bas_lib->basis_list) / sizeof(basis_t *);

    for (int i = 0; i < sz; i++) {
        basis_t *bas = bas_lib->basis_list[i];
        if (bas != NULL) {
            basis_delete(bas);
        }
    }

    free(bas_lib);
}


basis_t *basis_lib_set(basis_lib_t *bas_lib, basis_t *bas)
{
    int z = bas->element;
    if (z < 0 || z > N_CHEM_ELEMENTS) {
        return NULL;
    }

    bas_lib->basis_list[z] = bas;
    return bas;
}


basis_t *basis_lib_get(basis_lib_t *bas_lib, char *elem_symbol)
{
    int z = get_element_nuc_charge(elem_symbol);
    if (z == -1) {
        return NULL;
    }

    return bas_lib->basis_list[z];
}


void basis_lib_print(basis_lib_t *bas_lib)
{
    for (int i = 0; i < N_CHEM_ELEMENTS; i++) {
        if (bas_lib->basis_list[i] != NULL) {
            char buf[8];
            get_element_symbol(i, buf);
            printf(" ===================== %s ======================\n", buf);
            basis_print(bas_lib->basis_list[i]);
        }
    }
}
