#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ecp.h"


ecp_expansion_t *ecp_expansion_new(int nprim, int *powers, double *e, double *c)
{
    ecp_expansion_t *f = (ecp_expansion_t *) malloc(sizeof(ecp_expansion_t));

    f->powers = (int *) malloc(sizeof(int) * nprim);
    f->e = (double *) malloc(sizeof(double) * nprim);
    f->c = (double *) malloc(sizeof(double) * nprim);

    f->nprim = nprim;
    memcpy(f->powers, powers, sizeof(int) * nprim);
    memcpy(f->e, e, sizeof(double) * nprim);
    memcpy(f->c, c, sizeof(double) * nprim);

    return f;
}


void ecp_expansion_delete(ecp_expansion_t *f)
{
    free(f->powers);
    free(f->e);
    free(f->c);
    free(f);
}


ecp_t *ecp_new(int element)
{
    ecp_t *ecp = (ecp_t *) malloc(sizeof(ecp_t));

    ecp->element = element;
    ecp->n_core_elec = 0;
    memset(ecp->arep, 0, sizeof(ecp->arep));
    memset(ecp->esop, 0, sizeof(ecp->esop));

    return ecp;
}


void ecp_delete(ecp_t *ecp)
{
    int nfun = sizeof(ecp->arep) / sizeof(ecp_expansion_t *);

    for (int i = 0; i < nfun; i++) {
        if (ecp->arep[i]) {
            ecp_expansion_delete(ecp->arep[i]);
        }
        if (ecp->esop[i]) {
            ecp_expansion_delete(ecp->esop[i]);
        }
    }

    free(ecp);
}


void ecp_add_function(ecp_t *ecp, int type, int ang_mom, int nprim, int *powers, double *e, double *c)
{
    assert(type == ECP_AVERAGED || type == ECP_SPIN_ORBIT);

    ecp_expansion_t *f = ecp_expansion_new(nprim, powers, e, c);
    ecp_expansion_t **fun_list = (type == ECP_AVERAGED) ? ecp->arep : ecp->esop;

    assert(ang_mom <= ECP_MAX_ANG_MOM);

    int pos;
    if (ang_mom == ECP_UL) {
        pos = 0;
    }
    else {
        pos = ang_mom + 1;
    }

    if (fun_list[pos] != NULL) {
        ecp_expansion_delete(fun_list[pos]);
    }
    fun_list[pos] = f;
}


void ecp_print(ecp_t *ecp)
{
    printf("\nECP for E%d\n", ecp->element);
    printf(" %d core electrons\n", ecp->n_core_elec);
    printf(" arep:\n");
    for (int i = 0; i <= ECP_MAX_ANG_MOM; i++) {
        ecp_expansion_t *f = ecp->arep[i];
        if (f == NULL) {
            continue;
        }
        if (i == 0) {
            printf("L = Ul\n");
        }
        else {
            printf("L = %d\n", i - 1);
        }
        for (int j = 0; j < f->nprim; j++) {
            printf("%4d%18.8e%18.8e\n", f->powers[j], f->e[j], f->c[j]);
        }
    }
    printf(" esop:\n");
    for (int i = 0; i <= ECP_MAX_ANG_MOM; i++) {
        ecp_expansion_t *f = ecp->esop[i];
        if (f == NULL) {
            continue;
        }
        if (i == 0) {
            printf("L = Ul\n");
        }
        else {
            printf("L = %d\n", i - 1);
        }
        for (int j = 0; j < f->nprim; j++) {
            printf("%4d%18.8e%18.8e\n", f->powers[j], f->e[j], f->c[j]);
        }
    }
    printf("\n");
}
