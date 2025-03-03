#ifndef METRINIT_H
#define METRINIT_H

#include "galias.h"
#include <stdio.h>

void metr_init (int argc, const char *argv[]);
void print_helpc(FILE *out);
long process_argsc(int argc, const char *argv[], pkwd_sw &kwdsw);
long process_optargsc(int argc, const char *argv[], pkwd_sw &kwdsw);

#endif
