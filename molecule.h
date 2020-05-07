#ifndef MOLECULE_H_INCLUDED
#define MOLECULE_H_INCLUDED

#include "elements.h"

#define MAX_N_ATOMS 100

enum {
    SYMMETRY_AUTO,
    SYMMETRY_C1,
    SYMMETRY_Ci,
    SYMMETRY_Cs,
    SYMMETRY_C2,
    SYMMETRY_C2v,
    SYMMETRY_C2h,
    SYMMETRY_D2,
    SYMMETRY_D2h,
    SYMMETRY_Cinfv,
    SYMMETRY_Dinfh
};

typedef struct {
    int group;   // group label (one of SYMMETRY_* constants)
    int xyz[3];  // main group element (axis for axial groups and C2h, plane for Cs)
} symgrp_t;

typedef struct {
    int n_atoms;
    int charges[MAX_N_ATOMS];
    double x[MAX_N_ATOMS];
    double y[MAX_N_ATOMS];
    double z[MAX_N_ATOMS];
    symgrp_t sym_group;
} molecule_t;

molecule_t *molecule_new();
void molecule_delete(molecule_t *mol);
void molecule_add_atom(molecule_t *mol, int nuc_charge, double x, double y, double z);
void molecule_print(molecule_t *mol);
int molecule_n_atom_types(molecule_t *mol);
int molecule_n_atoms_of(molecule_t *mol, int element);

#endif /* MOLECULE_H_INCLUDED */
