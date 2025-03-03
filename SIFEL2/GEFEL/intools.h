#ifndef INTOOLS_H
#define INTOOLS_H

#include <stdio.h>
#include "xfile.h"

long isgkw(char *s, char *&ptr);
long ispkw(char *s, char *&ptr);
long iscomment(char *s, char *&ptr);
void delcomment(char *s);
FILE *cleancommblock(FILE *in);
void geteoln(FILE *in);
void getnextln(FILE *in);
void getnexttxt(FILE *in);
long getgkwid(FILE *in, long id);
long getgkwid(XFILE *in, long id, long &numl);
long getpkw(FILE *in);
long getnprop(FILE *in);
long getnprop(XFILE *in, long &numl);
long getlong(FILE *in, long &num);
long getint(FILE *in, int &num);
long getdouble(FILE *in, double &num);
long getstring(FILE *in, char *str, long maxlen);
long getstring2(FILE *in, char *str, long maxlen);
long inputln(FILE* in, char* anystr, long maxlen);

#endif
