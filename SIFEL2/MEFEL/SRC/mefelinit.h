#ifndef MEFELINIT_H
#define MEFELINIT_H

#include <stdio.h>
#include "galias.h"

/// function initializes 
void mefel_init (int argc, const char *argv[]);
/// function prints help message to the stdout
void print_help(FILE *out, const char *prgname);
/// function decodes used command line arguments
long process_args(int argc, const char *argv[], pkwd_sw &kwdsw);

#endif
