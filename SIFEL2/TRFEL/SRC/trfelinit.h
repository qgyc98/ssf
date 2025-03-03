#ifndef TRFELINIT_H
#define TRFELINIT_H

#include "galias.h"

void trfel_init (int argc, const char *argv[]);
void print_helpt(FILE *out, const char *prgname);
long process_argst(int argc, const char *argv[], pkwd_sw &kwdsw);

void radiation_init ();

#endif
