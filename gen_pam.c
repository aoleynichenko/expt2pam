#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basis_lib.h"
#include "ecp_lib.h"
#include "elements.h"
#include "molecule.h"

#define MAX_ANG_MOM 10

#define MAX_CONTRACTED 100

void basis_get_blocks_for_dirac(basis_t *bas, int *n_blocks, int* n_block_sizes, int len_block_sizes);
char angular_momentum_to_char(int l);
int basis_max_ang_mom(basis_t *bas);
void print_l_blocks_dirac(FILE *out, basis_t *bas, int L);

int basis_get_number_l_blocks(basis_t *bas, int L);

void gen_pam(FILE *out, molecule_t *mol, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib)
{
    // first three lines -- comment
    fprintf(out, "\n");
    fprintf(out, "\n");
    fprintf(out, "\n");

    // 'C' for Cartesian GTO basis;
    // number of atom types
    // 'A' label for angstroms if needed
    int n_atom_types = molecule_n_atom_types(mol);
    fprintf(out, "C  %2d\n", n_atom_types);

    // loop over atom types
    // heavier elements first
    for (int i = N_CHEM_ELEMENTS-1; i >= 0; i--) {
        char elem_sym[MAX_ELEMENT_SYMBOL];

        int n_atoms = molecule_n_atoms_of(mol, i);
        if (n_atoms == 0) {
            continue;
        }

        get_element_symbol(i, elem_sym);

        // generate XYZ coordinates for atoms of the given type
        fprintf(out, "      %3d.   %2d\n", i, n_atoms);
        for (int j = 0, k = 1; j < mol->n_atoms; j++) {
            if (mol->charges[j] == i) {
                fprintf(out, "%s%d%16.10f%16.10f%16.10f\n", elem_sym, k++, mol->x[j], mol->y[j], mol->z[j]);
            }
        }

        // generate basis
        basis_t *bas = basis_lib_get(bas_lib, elem_sym);
        int n_blocks;
        int block_sizes[MAX_ANG_MOM];

        int max_ang_mom = basis_max_ang_mom(bas);

        fprintf(out, "LARGE EXPLICIT%5d", max_ang_mom+1);
        for (int l = 0; l <= max_ang_mom; l++) {
            int n_blocks_dirac = basis_get_number_l_blocks(bas, l);
            fprintf(out, "%5d", n_blocks_dirac);
        }
        fprintf(out, "\n");

        for (int l = 0; l <= max_ang_mom; l++) {
            print_l_blocks_dirac(stdout, bas, l);
        }

        // print ECP if needed
        ecp_t *ecp = ecp_lib_get(ecp_lib, elem_sym);
        if (ecp == NULL) {
            continue;
        }
        int len_arep, len_esop;
        ecp_get_len(ecp, &len_arep, &len_esop);
        fprintf(out, "ECP %d %d %d\n", ecp->n_core_elec, len_arep, len_esop);

        for (int i = 0; i <= ECP_MAX_ANG_MOM; i++) {
            ecp_expansion_t *f = ecp->arep[i];
            if (f == NULL) {
                continue;
            }
            if (i == 0) {
                fprintf(out, "# ul\n");
            }
            else {
                fprintf(out, "# %c\n", tolower(angular_momentum_to_char(i - 1)));
            }
            for (int j = 0; j < f->nprim; j++) {
                fprintf(out, "%4d%24.12e%24.12e\n", f->powers[j], f->e[j], f->c[j]);
            }
        }
        for (int i = 0; i <= ECP_MAX_ANG_MOM; i++) {
            ecp_expansion_t *f = ecp->esop[i];
            if (f == NULL) {
                continue;
            }
            else {
                fprintf(out, "# %c (spin-orbit)\n", tolower(angular_momentum_to_char(i - 1)));
            }
            for (int j = 0; j < f->nprim; j++) {
                fprintf(out, "%4d%24.12e%24.12e\n", f->powers[j], f->e[j], f->c[j]);
            }
        }
    }

    fprintf(out, "FINISH\n");
}


int basis_max_ang_mom(basis_t *bas)
{
    int max_ang_mom = -1;

    for (int ibf = 0; ibf < bas->nfun; ibf++) {
        bfn_t *bf = bas->functions[ibf];
        if (bf->l > max_ang_mom) {
            max_ang_mom = bf->l;
        }
    }

    return max_ang_mom;
}


/**
 * returns number of L-blocks for given angular momentum L
 */
int basis_get_number_l_blocks(basis_t *bas, int L)
{
    bfn_t *first_basis_functions[MAX_CONTRACTED];   // functions which are first in each block
    int block_sizes[MAX_CONTRACTED];
    int n_blocks = 0;
    int has_primitive_functions = 0;

    memset(first_basis_functions, 0, sizeof(first_basis_functions));
    memset(block_sizes, 0, sizeof(block_sizes));

    for (int ibf = 0; ibf < bas->nfun; ibf++) {
        bfn_t *bf = bas->functions[ibf];
        if (bf->l != L) {
            continue;
        }

        // find block for this function
        int ib;
        if (bf->nprim == 1) {
            has_primitive_functions = 1;
            continue;
        }
        for (ib = 0; ib < n_blocks; ib++) {
            if (bfn_same_exponents(first_basis_functions[ib], bf)) {
                block_sizes[ib]++;
                break;
            }
        }
        if (ib == n_blocks) { // block for this function was not found
            first_basis_functions[n_blocks] = bf;
            block_sizes[n_blocks] = 1;
            n_blocks++;
        }
    }

    int n_blocks_dirac = 0;
    for (int i = 0; i < n_blocks; i++) {
        n_blocks_dirac += block_sizes[i] / 4 + (block_sizes[i] % 4 != 0);
    }
    n_blocks_dirac += has_primitive_functions;

    return n_blocks_dirac;
}


void print_l_blocks_dirac(FILE *out, basis_t *bas, int L)
{
    bfn_t *basis_functions[MAX_CONTRACTED][MAX_CONTRACTED];   // functions which are first in each block
    int block_sizes[MAX_CONTRACTED];
    int n_blocks = 0;
    int has_primitive_functions = 0;
    int n_primitive_functions = 0;

    memset(basis_functions, 0, sizeof(basis_functions));
    memset(block_sizes, 0, sizeof(block_sizes));

    for (int ibf = 0; ibf < bas->nfun; ibf++) {
        bfn_t *bf = bas->functions[ibf];
        if (bf->l != L) {
            continue;
        }

        // find block for this function
        int ib;
        if (bf->nprim == 1) {
            has_primitive_functions = 1;
            n_primitive_functions++;
            continue;
        }
        for (ib = 0; ib < n_blocks; ib++) {
            if (bfn_same_exponents(basis_functions[ib][0], bf)) { // add function to the existing block
                int bs = block_sizes[ib];
                basis_functions[ib][bs] = bf;
                block_sizes[ib]++;
                break;
            }
        }
        if (ib == n_blocks) { // block for this function was not found; create new block
            basis_functions[n_blocks][0] = bf;
            block_sizes[n_blocks] = 1;
            n_blocks++;
        }
    }

    if (n_blocks == 0 && has_primitive_functions == 0) {
        return;
        // nothing to print
    }

    fprintf(out, "# %c functions\n", tolower(angular_momentum_to_char(L)));

    for (int ibl = 0; ibl < n_blocks; ibl++) {
        int block_size = block_sizes[ibl];
        bfn_t **funcs = basis_functions[ibl];
        int nprim = funcs[0]->nprim;
        fprintf(out, "F%3d%3d\n", nprim, block_size);
        for (int i = 0; i < nprim; i++) {
            fprintf(out, "%24.10f", funcs[0]->e[i]);
            for (int j = 0; j < block_size; j++) {
                fprintf(out, "%16.10f", funcs[j]->c[i]);
            }
            fprintf(out, "\n");
        }
    }

    // print primitive Gaussians
    if (has_primitive_functions) {
        fprintf(out, "F%3d%3d\n", n_primitive_functions, 0);
        for (int ibf = 0; ibf < bas->nfun; ibf++) {
            bfn_t *bf = bas->functions[ibf];
            if (bf->l != L || bf->nprim != 1) {
                continue;
            }
            fprintf(out, "%24.10f\n", bf->e[0]);
        }
    }
}


void basis_get_blocks_for_dirac(basis_t *bas, int *n_blocks, int* n_block_sizes, int len_block_sizes)
{
    int n_cont_fun_by_l[MAX_ANG_MOM];
    int n_prim_fun_by_l[MAX_ANG_MOM];

    memset(n_cont_fun_by_l, 0, sizeof(n_cont_fun_by_l));
    memset(n_prim_fun_by_l, 0, sizeof(n_prim_fun_by_l));

    for (int i = 0; i < len_block_sizes; i++) {
        n_block_sizes[i] = 0;
    }

    // count functions in each irrep of SO(3) (L)
    for (int ifun = 0; ifun < bas->nfun; ifun++) {
        bfn_t *bf = bas->functions[ifun];
        int ang_mom = bf->l;
        if (bf->nprim == 1) {
            n_prim_fun_by_l[ang_mom]++;
        }
        else {
            n_cont_fun_by_l[ang_mom]++;
        }
    }

    // for each angulam momentum: max 4 contracted functions in block
    // + 1 separate block for uncontracted Gaussians
    for (int i = 0; i < len_block_sizes; i++) {
        n_block_sizes[i] = n_cont_fun_by_l[i] / 4 + (n_cont_fun_by_l[i] % 4 != 0) + (n_prim_fun_by_l[i] > 0);
    }

    *n_blocks = 0;
    for (int i = 0; i < len_block_sizes; i++) {
        if (n_block_sizes[i] != 0) {
            *n_blocks = i+1;
        }
    }
}





char angular_momentum_to_char(int l)
{
    char labels[] = "SPDFGHIKLM";

    assert(l >=0 && l < (sizeof(labels) / sizeof(char)));

    return labels[l];
}
