#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "molecule.h"

molecule_t *molecule_new()
{
    molecule_t *mol = (molecule_t *) malloc(sizeof(molecule_t));

    mol->n_atoms = 0;
    mol->sym_group.group = SYMMETRY_AUTO;
    mol->sym_group.xyz[0] = 0;
    mol->sym_group.xyz[1] = 0;
    mol->sym_group.xyz[2] = 0;

    return mol;
}


void molecule_delete(molecule_t *mol)
{
    free(mol);
}


void molecule_add_atom(molecule_t *mol, int nuc_charge, double x, double y, double z)
{
    int n_at = mol->n_atoms;

    assert(nuc_charge >= 0 && nuc_charge < N_CHEM_ELEMENTS);

    mol->charges[n_at] = nuc_charge;
    mol->x[n_at] = x;
    mol->y[n_at] = y;
    mol->z[n_at] = z;
    mol->n_atoms++;
}


int molecule_n_atom_types(molecule_t *mol)
{
    int n_atoms[N_CHEM_ELEMENTS];
    int n_atom_types;

    memset(n_atoms, 0, sizeof(n_atoms));

    for (int i = 0; i < mol->n_atoms; i++) {
        n_atoms[mol->charges[i]]++;
    }

    n_atom_types = 0;
    for (int i = 0; i < N_CHEM_ELEMENTS; i++) {
        if (n_atoms[i] > 0) {
            n_atom_types++;
        }
    }

    return n_atom_types;
}


int molecule_n_atoms_of(molecule_t *mol, int element)
{
    int n_atoms = 0;

    for (int i = 0; i < mol->n_atoms; i++) {
        if (mol->charges[i] == element) {
            n_atoms++;
        }
    }

    return n_atoms;
}


void molecule_print(molecule_t *mol)
{
    printf("\ngeometry:\n");
    printf("---------\n");
    for (int i = 0; i < mol->n_atoms; i++) {
        printf("  %3d%12.8f%12.8f%12.8f\n", mol->charges[i], mol->x[i], mol->y[i], mol->z[i]);
    }
    printf("Symmetry: %d\n (orientation %d %d %d)\n", mol->sym_group.group,
        mol->sym_group.xyz[0], mol->sym_group.xyz[1], mol->sym_group.xyz[2]);
    printf("\n");
}
