#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basis_lib.h"
#include "ecp_lib.h"
#include "expt.h"
#include "molecule.h"

#define MAX_FILE_NAME 1024

void parse_argv(int argc, char **argv, char *file_name, int maxlen);


int main(int argc, char **argv)
{
    char inp_file[MAX_FILE_NAME];
    basis_lib_t *bas;
    ecp_lib_t *ecp;
    molecule_t *mol;

    parse_argv(argc, argv, inp_file, MAX_FILE_NAME);

    bas = basis_lib_new();
    ecp = ecp_lib_new();
    mol = molecule_new();

    expt_parse(inp_file, mol, bas, ecp);

    /*basis_lib_print(bas);
    ecp_lib_print(ecp);
    molecule_print(mol);*/

    gen_pam(stdout, mol, bas, ecp);

    molecule_delete(mol);
    ecp_lib_delete(ecp);
    basis_lib_delete(bas);
}


void parse_argv(int argc, char **argv, char *file_name, int maxlen)
{
    if (argc != 2) {
        printf("Usage: %s <input-file>\n", argv[0]);
        exit(1);
    }

    strncpy(file_name, argv[1], maxlen);
    file_name[maxlen-1] = '\0';
}
