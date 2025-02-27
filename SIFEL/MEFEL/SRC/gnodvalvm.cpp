#include "gnodvalvm.h"

#include <stdlib.h>

gnodvalvm::gnodvalvm() : displv{NULL}, loadv{NULL}, iforv{NULL}, residv{NULL}
{
}



gnodvalvm::gnodvalvm(double *d, double *fl, double *fi, double *fr) : displv{d}, loadv{fl}, iforv{fi}, residv{fr}
{
}
