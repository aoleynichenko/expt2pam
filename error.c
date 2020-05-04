#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "error.h"

#define MAX_ERR_LEN 1024


// writes error message & aborts execution
void errquit(char *fmt, ...)
{
    va_list args;
    char errbuf[MAX_ERR_LEN];

    va_start(args, fmt);
    vsprintf(errbuf, fmt, args);
    errbuf[MAX_ERR_LEN - 1] = '\0';
    va_end(args);
    printf("\nERROR:\n");
    printf("%s\n", errbuf);

    // abort execution
    exit(1);
}
