#ifndef EXPT_H_INCLUDED
#define EXPT_H_INCLUDED

#include "basis.h"
#include "basis_lib.h"
#include "molecule.h"

void expt_parse(char *file_name, molecule_t *mol, basis_lib_t *basis, ecp_lib_t *ecp_lib);
void gen_pam(FILE *out, molecule_t *mol, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib);

#endif /* EXPT_H_INCLUDED */
