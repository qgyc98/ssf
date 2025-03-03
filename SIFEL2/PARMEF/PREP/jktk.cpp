#include "../SRC/pglobal.h"
#include "../SRC/seqfile.h"
#include "../SRC/paral.h"
#include "../../MEFEL/PREP/input.h"
#include <stdio.h>
#include <string.h>
#include "mpi.h"

probdesc *Mp;
pprobdesc *Pmp;
mechtop *Mt;
paral *Plg;
mechmat *Mm;
mechcrsec *Mc;
mechbclc *Mb;

densemat *Sm_dm,*Mm_dm;
skyline *Sm_sky,*Mm_sky;
comprow *Sm_cr,*Mm_cr;
symcomprow *Sm_scr,*Mm_scr;
lhsrhs *Lsrs;

barel2d *Bar2d=NULL;
planestresslt *Pslt=NULL;
planestressqt *Psqt=NULL;
planestressrotlt *Psrlt=NULL;
planestresslq *Pslq=NULL;
planestressqq *Psqq=NULL;
planestressrotlq *Psrlq=NULL;
cctelem *Cct=NULL;
axisymlq *Asymlq=NULL;
lintet *Ltet=NULL;
linhex *Lhex=NULL;
quadhex *Qhex=NULL;

int Nproc,Myrank,Ndom;
long Ndof,Indof,Mespr;
long *Domproc;
FILE *Out;

int main (int argc,char *argv[])
{
  FILE *in, *tmpf;
  descrip *d;

  fprintf (stderr,"\n\n *** PARMEFEL PREPROCESSOR ***\n");
  fprintf (stderr," -----------------------------\n");

  Mp  = new probdesc;
  Pmp = new pprobdesc;
  Mt  = new mechtop;
  Mm  = new mechmat;
  Mc  = new mechcrsec;
  Mb  = new mechbclc;

  if (argc < 3){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : mechprep input_file_name output_file name\n\n");
    delete Mp;  delete Mt;  delete Mm;  delete Mc;  delete Mb;
    return(1);
  }
  in = fopen(argv[1],"r");
  if (in==NULL){
    fprintf (stderr,"\n Input file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    delete Mp;  delete Mt;  delete Mm;  delete Mc;  delete Mb;
    return(2);
  }

  d = new descrip;
  // reading line with topology file name
  getstring(in, d->t3df, 1024);
  // reading line with material database file name
  getstring(in, d->matf, 1024);
  // reading line with cross-section database file name
  getstring(in, d->crf, 1024);
  // reading line with topology file format indicator,
  getlong(in,d->t3d);
  d->seq = 0;
  // reading problem description
  tmpf = cleancommblock(in);
  if (tmpf == NULL)
  {
    fprintf(stderr, "\n\nError - unable to open temporary file");
    fprintf(stderr, "\n in function main() (mechprep.cpp)\n");
    return (1);
  }
  Pmp->read(tmpf);
  fclose(tmpf);
  Mp->tprob = Pmp->tprob;
  Mp->nlcsip = Pmp->nlcsip;
  Mp->nstrcomp = Pmp->nstrcomp;

  input(in, d);
  fclose(in);

  FILE *out = fopen(argv[2], "wt");
  if (out==NULL){
    fprintf (stderr,"\n Output file couldn't be opened.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(7);
  }
  Pmp->print(out);
  output(out, d);
  fclose(out);

  fprintf(stderr, "\n");

  delete d;
  delete Pmp;  delete Mp;  delete Mt;  delete Mm;  delete Mc;  delete Mb;
  return(0);
}
