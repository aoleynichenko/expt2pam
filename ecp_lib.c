#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ecp.h"
#include "ecp_lib.h"
#include "error.h"

ecp_lib_t *ecp_lib_new()
{
    ecp_lib_t *ecp_lib = (ecp_lib_t *) malloc(sizeof(ecp_lib_t));

    memset(ecp_lib->ecp_list, 0, sizeof(ecp_lib->ecp_list));

    return ecp_lib;
}


void ecp_lib_delete(ecp_lib_t *ecp_lib)
{
    int sz = sizeof(ecp_lib->ecp_list) / sizeof(ecp_t *);

    for (int i = 0; i < sz; i++) {
        ecp_t *ecp = ecp_lib->ecp_list[i];
        if (ecp != NULL) {
            ecp_delete(ecp);
        }
    }

    free(ecp_lib);
}


ecp_t *ecp_lib_set(ecp_lib_t *ecp_lib, ecp_t *ecp)
{
    int z = ecp->element;
    if (z < 0 || z > N_CHEM_ELEMENTS) {
        return NULL;
    }

    ecp_lib->ecp_list[z] = ecp;
    return ecp;
}


ecp_t *ecp_lib_get(ecp_lib_t *ecp_lib, char *elem_symbol)
{
    int z = get_element_nuc_charge(elem_symbol);
    if (z == -1) {
        return NULL;
    }

    return ecp_lib->ecp_list[z];
}


void ecp_lib_print(ecp_lib_t *ecp_lib)
{
    for (int i = 0; i < N_CHEM_ELEMENTS; i++) {
        if (ecp_lib->ecp_list[i] != NULL) {
            char buf[8];
            get_element_symbol(i, buf);
            printf(" ===================== %s ======================\n", buf);
            ecp_print(ecp_lib->ecp_list[i]);
        }
    }
}
