#ifndef ECP_H_INCLUDED
#define ECP_H_INCLUDED

#include "elements.h"

#define ECP_MAX_ANG_MOM 10

#define ECP_UL -1

#define ECP_AVERAGED    0
#define ECP_SPIN_ORBIT  1

typedef struct {
    int nprim;
    int *powers;
    double *e;
    double *c;
} ecp_expansion_t;

typedef struct {
    int element;
    int n_core_elec;
    ecp_expansion_t *arep[ECP_MAX_ANG_MOM+1]; // Ul in [0], S in [1], etc
    ecp_expansion_t *esop[ECP_MAX_ANG_MOM+1];
} ecp_t;


ecp_expansion_t *ecp_expansion_new(int nprim, int *powers, double *e, double *c);
void ecp_expansion_delete(ecp_expansion_t *ecp);
ecp_t *ecp_new(int element);
void ecp_delete(ecp_t *ecp);
void ecp_add_function(ecp_t *ecp, int type, int ang_mom, int nprim, int *powers, double *e, double *c);
void ecp_print(ecp_t *ecp);
void ecp_get_len(ecp_t *ecp, int *len_arep, int *len_esop);


typedef struct {
    ecp_t *ecp_list[N_CHEM_ELEMENTS];
} ecp_lib_t;

ecp_lib_t *ecp_lib_new();
void ecp_lib_delete(ecp_lib_t *);
ecp_t *ecp_lib_set(ecp_lib_t *ecp_lib, ecp_t *ecp);
ecp_t *ecp_lib_get(ecp_lib_t *ecp_lib, char *elem_symbol);
void ecp_lib_print(ecp_lib_t *ecp_lib);


#endif /* ECP_H_INCLUDED */
