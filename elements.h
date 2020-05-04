#ifndef ELEMENTS_H_INCLUDED
#define ELEMENTS_H_INCLUDED

#define N_CHEM_ELEMENTS 131 // 130 + Ghost
#define MAX_ELEMENT_SYMBOL 8

typedef struct {
    int z;
    char sym[MAX_ELEMENT_SYMBOL];
} element_t;

void get_element_symbol(int z, char *sym);
int get_element_nuc_charge(char *sym);

#endif /* ELEMENTS_H_INCLUDED */
