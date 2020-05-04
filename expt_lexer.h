#ifndef EXPT_LEXER_H_INCLUDED
#define EXPT_LEXER_H_INCLUDED

#include <stdio.h>

extern FILE *yyin;
extern char *yytext;
extern int yylineno;
extern int yycol;

int yylex();

enum {
    /* token types */
    TT_QUOTE = 0,
    TT_STAR,
    TT_HYPHEN,
    TT_INTEGER,
    TT_FLOAT,
    TT_WORD,
    /* keywords */
    KEYWORD_GEOMETRY,
    KEYWORD_BASIS,
    KEYWORD_ECP,
    KEYWORD_END,
    /* special */
    END_OF_LINE,
    END_OF_FILE
};

#endif /* EXPT_LEXER_H_INCLUDED */
