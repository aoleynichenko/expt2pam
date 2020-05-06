#ifndef ECP_LIB_H_INCLUDED
#define ECP_LIB_H_INCLUDED

#include "ecp.h"
#include "elements.h"

typedef struct {
    ecp_t *ecp_list[N_CHEM_ELEMENTS];
} ecp_lib_t;

ecp_lib_t *ecp_lib_new();
void ecp_lib_delete(ecp_lib_t *);
ecp_t *ecp_lib_set(ecp_lib_t *ecp_lib, ecp_t *ecp);
ecp_t *ecp_lib_get(ecp_lib_t *ecp_lib, char *elem_symbol);
void ecp_lib_print(ecp_lib_t *ecp_lib);

#endif /* ECP_LIB_H_INCLUDED */
