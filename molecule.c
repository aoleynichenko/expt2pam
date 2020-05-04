#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "molecule.h"

molecule_t *molecule_new()
{
    molecule_t *mol = (molecule_t *) malloc(sizeof(molecule_t));

    mol->n_atoms = 0;

    return mol;
}


void molecule_delete(molecule_t *mol)
{

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


void molecule_print(molecule_t *mol)
{
    printf("\ngeometry:\n");
    printf("---------\n");
    for (int i = 0; i < mol->n_atoms; i++) {
        printf("  %3d%12.8f%12.8f%12.8f\n", mol->charges[i], mol->x[i], mol->y[i], mol->z[i]);
    }
    printf("\n");
}
