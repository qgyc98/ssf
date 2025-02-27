#ifndef STACKTRACE_H
#define STACKTRACE_H
#include <stdio.h>

/// function sets name of the program
void set_prgname(const char *name);
/// function prints stack trace to the file out
void stack_trace(FILE *out, long level);
/// function for stack trace in Linux
void lin_stack_trace(FILE *out, long level); 
/// function performs name demangling for Watcom compiler
void demangle_watcom(const char *mname, char *dname);
/// function reduces path length of the source file name
char *remove_path_level(char *path, long level);
/// function for stack trace in Windows
void win_stack_trace(FILE *out, long level); 
#endif
