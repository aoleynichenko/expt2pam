#ifndef MOLECULE_H_INCLUDED
#define MOLECULE_H_INCLUDED

#include "elements.h"

#define MAX_N_ATOMS 100

typedef struct {
    int n_atoms;
    int charges[MAX_N_ATOMS];
    double x[MAX_N_ATOMS];
    double y[MAX_N_ATOMS];
    double z[MAX_N_ATOMS];
} molecule_t;

molecule_t *molecule_new();
void molecule_delete(molecule_t *mol);
void molecule_add_atom(molecule_t *mol, int nuc_charge, double x, double y, double z);
void molecule_print(molecule_t *mol);

#endif /* MOLECULE_H_INCLUDED */
