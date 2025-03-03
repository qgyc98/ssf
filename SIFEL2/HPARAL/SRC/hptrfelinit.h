#ifndef HPTRFELINIT_H
#define HPTRFELINIT_H

struct XFILE;

void hptrfel_init (int argc,const char *argv[]);

void hptrfel_init_tiles (int argc,const char *argv[]);

void hptrfel_init_tiles_parsolver (int argc,const char *argv[], XFILE* &in);

#endif
