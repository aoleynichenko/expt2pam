#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basis.h"
#include "basis_lib.h"
#include "error.h"
#include "expt_lexer.h"
#include "molecule.h"

#define MAX_LINE_LEN 1024
#define MAX_FILE_NAME 1024
#define MAX_PRIMITIVES 100
#define MAX_CONTRACTED 100

void directive_geometry(molecule_t *mol);
void directive_basis(basis_lib_t *basis);
void directive_basis_block(basis_lib_t *basis);

void yyerror(char *s);
int next_token();
void put_back(int token_type);
int match(int required_type);
void str_tolower(char *s);
int parse_state_spec(char *buf, char *rep_name, int *state_no);
void set_input_file_name(char *file_name);
void parse_sector(char *s, int *h, int *p);
int angular_momentum_to_int(char *lstr);


void expt_parse(char *file_name, molecule_t *mol, basis_lib_t *basis)
{
    int token_type;

    yyin = fopen(file_name, "r");
    if (yyin == NULL) {
        errquit("Input file \"%s\" not found", file_name);
    }

    set_input_file_name(file_name);

    token_type = next_token();
    while (token_type != END_OF_FILE) {
        switch (token_type) {
            case KEYWORD_GEOMETRY:
                directive_geometry(mol);
                break;
            case KEYWORD_BASIS:
                directive_basis(basis);
                break;
            case END_OF_LINE:
                // nothing to do
                break;
            default:
                yyerror("unknown keyword");
                break;
        }

        /* each "non-newline" directive must end with EOL or EOF */
        if (token_type != END_OF_LINE && token_type != END_OF_FILE) {
            token_type = next_token();
            if (token_type != END_OF_LINE && token_type != END_OF_FILE) {
                yyerror("end of line is expected");
            }
        }

        /* go to the next directive */
        token_type = next_token();
    }

    fclose(yyin);
}


void directive_geometry(molecule_t *mol)
{
    static char *err_no_newline = "expected end of line";
    static char *err_no_elem    = "wrong element symbol";
    static char *err_no_float   = "float number is expected";
    int token_type;

    token_type = next_token();
    if (token_type != END_OF_LINE) {
        yyerror(err_no_newline);
    }

    while (1) {
        token_type = next_token();

        if (token_type == TT_WORD) { // add new atom
            int nuc_charge;
            double r[] = {0.0, 0.0, 0.0};

            // read element symbol -> nuclear charge
            nuc_charge = get_element_nuc_charge(yytext);
            if (nuc_charge == -1) {
                yyerror(err_no_elem);
            }

            // read nuclear coordinates x,y,z
            for (int i = 0; i < 3; i++) {
                token_type = next_token();
                if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
                    yyerror(err_no_float);
                }
                r[i] = atof(yytext);
            }

            token_type = next_token();
            if (token_type != END_OF_LINE) {
                yyerror(err_no_newline);
            }

            molecule_add_atom(mol, nuc_charge, r[0], r[1], r[2]);
        }
        else if (token_type == END_OF_LINE) {
            /* empty line, nothing to do */
            continue;
        }
        else if (token_type == KEYWORD_END) {
            return;
        }
        else {
            yyerror("unexpected token");
        }
    }

}


void directive_basis(basis_lib_t *basis)
{
    static char *err_no_newline = "expected end of line";
    int token_type;

    token_type = next_token();
    if (token_type != END_OF_LINE) {
        yyerror(err_no_newline);
    }

    // read blocks of functions until the 'end' keyword

    while (1) {
        token_type = next_token();

        if (token_type == TT_WORD) {
            put_back(token_type);
            directive_basis_block(basis);
        }
        else if (token_type == END_OF_LINE) {
            /* empty line, nothing to do */
            continue;
        }
        else if (token_type == KEYWORD_END) {
            return;
        }
        else {
            yyerror("unexpected token");
        }
    }
}


void directive_basis_block(basis_lib_t *basis_lib)
{
    static char *err_ang_mom   = "angular momentum symbol (one of SPDFGHIKLM) is expected";
    static char *err_end_line  = "end of line is expected";
    static char *err_end_float = "float number or end of line are expected";
    static char *err_float     = "float number if expected";
    static char *err_no_fun    = "no functions found in this block";
    static char *err_no_elem   = "wrong element symbol";

    int token_type;
    static double matrix[MAX_CONTRACTED][MAX_PRIMITIVES];
    int n_primitives = 0;
    int n_contracted = 0;
    basis_t *basis = NULL;

    memset(matrix, 0, sizeof(matrix));

    token_type = next_token();
    str_tolower(yytext);
    basis = basis_lib_get(basis_lib, yytext);
    if (basis == NULL) {
        int z = get_element_nuc_charge(yytext);
        if (z == -1) {
            yyerror(err_no_elem);
        }
        basis = basis_new(z);
    }

    token_type = next_token();
    if (token_type != TT_WORD) {
        yyerror(err_ang_mom);
    }
    int ang_mom = angular_momentum_to_int(yytext);
    if (ang_mom == -1) {
        yyerror(err_ang_mom);
    }

    token_type = next_token();
    if (token_type != END_OF_LINE) {
        yyerror(err_end_line);
    }

    while (1) { // loop over rows
        token_type = next_token();
        if (token_type == END_OF_LINE) { // just empty line
            continue;
        }
        else if (token_type == TT_WORD || token_type == KEYWORD_END) { // end of block
            put_back(token_type);
            break;
        }
        else if (token_type == TT_INTEGER || token_type == TT_FLOAT) {
            // loop over columns
            double exponent = atof(yytext);
            matrix[0][n_primitives] = exponent;
            n_contracted = 0;
            while ((token_type = next_token()) == TT_INTEGER || token_type == TT_FLOAT) {
                double coef = atof(yytext);
                matrix[1+n_contracted][n_primitives] = coef;
                n_contracted++;
            }
            if (token_type != END_OF_LINE) {
                yyerror(err_end_float);
            }
            n_primitives++;
        }
        else {
            yyerror(err_float);
        }
    }

    if (n_primitives == 0) {
        yyerror(err_no_fun);
    }

    if (n_contracted == 0) { // special case of uncontracted functions
        for (int i = 0; i < n_primitives; i++) {
            double coef_1 = 1.0;
            basis_add_function(basis, ang_mom, 1, &matrix[0][i], &coef_1);
        }
    }
    else {
        for (int i = 0; i < n_contracted; i++) {
            basis_add_function(basis, ang_mom, n_primitives, matrix[0], matrix[1+i]);
        }
    }

    basis_lib_set(basis_lib, basis);
}


/**
 * Converts string containing angular momentum symbol (S,P,D,...) to its integer
 * value: S -> 0, P -> 1, ...
 * returns -1 in case of error.
 */
int angular_momentum_to_int(char *lstr)
{
    if (strlen(lstr) != 1) {
        return -1;
    }

    switch (lstr[0]) {
        case 'S': case 's':
            return 0;
        case 'P': case 'p':
            return 1;
        case 'D': case 'd':
            return 2;
        case 'F': case 'f':
            return 3;
        case 'G': case 'g':
            return 4;
        case 'H': case 'h':
            return 5;
        case 'I': case 'i':
            return 6;
        case 'K': case 'k':
            return 7;
        case 'L': case 'l':
            return 8;
        case 'M': case 'm':
            return 9;
        default:
            return -1;
    }
}


/*******************************************************************************
                     PARSER OF INPUT FILES - GLOBAL VARIABLES
 ******************************************************************************/

// flag: was the last token returned to the lexer?
static int pushed_back = 0;

// name of the input file with language to be processed
static char input_file_name[MAX_FILE_NAME];


/*******************************************************************************
                     PARSER OF INPUT FILES - HELPER FUNCTIONS
 ******************************************************************************/


void set_input_file_name(char *file_name)
{
    strcpy(input_file_name, file_name);
}


void print_line_by_num(char *file_name, int required_lineno)
{
    FILE *f;
    int c;
    int lineno;

    f = fopen(file_name, "r");
    if (f == NULL) {
        errquit("print_line_by_num(): cannot open file '%s'", file_name);
    }

    lineno = 1;
    while ((c = fgetc(f)) != EOF) {
        if (lineno == required_lineno) {
            printf(" ");
            do {
                putchar(c);
                c = fgetc(f);
            }
            while (c != '\n' && c != EOF);
            printf("\n");
            fclose(f);
            return;
        }
        if (c == '\n') {
            lineno++;
        }
    }

    errquit("print_line_by_num(): no line number %d in '%s'", required_lineno, file_name);
    fclose(f);
}

void yyerror(char *s)
{
    int lineno = yylineno;

    // lineno = yylineno-1 if last token was END_OF_LINE
    if (strcmp(yytext, "\n") == 0) {
        lineno = yylineno - 1;
    }
    fclose(yyin);

    printf("%s:%d:%d: error: %s\n\n", input_file_name, lineno, yycol, s);
    print_line_by_num(input_file_name, lineno);

    // print pointer to the erroneous token
    printf(" ");
    for (int i = 0; i < yycol; i++) {
        printf(" ");
    }
    printf("^\n\n");
    printf("execution terminated.\n");

    exit(1);
}


int next_token()
{
    if (pushed_back) {
        int token_type = pushed_back;
        pushed_back = 0;
        return token_type;
    }
    return yylex();
}


void put_back(int token_type)
{
    pushed_back = token_type;
}

int match(int required_type)
{
    int token_type = next_token();
    if (token_type != required_type) {
        return 0;
    }
    return 1;
}


void str_tolower(char *s)
{
    for ( ; *s; s++ ) {
        *s = tolower(*s);
    }
}
