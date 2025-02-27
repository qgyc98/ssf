#ifndef LIBTRACE_H
#define LIBTRACE_H

#include <string.h>

/// sets initializes internal variable with the program name
void set_program_name(const char *name);

/// returns pointer to string with the program name
const char *get_program_name();

/// returns pointer of bfd structure
void *get_bfd_pointer();

/// initializes system of detection of source files and lines
int libtrace_init(const char *file_name, const char *section_name, const char *target);

/// returns source file name and line (optionally also function name) for given code address
int libtrace_resolve(void *addr, char *buf_func, size_t buf_func_len, 
                     char *buf_file, size_t buf_file_len, unsigned int &line);

/// closes system of detection of source files and lines
void libtrace_close(void);

#endif
