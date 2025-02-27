#define EXTERN
#include "globalt.h"

#include "stochdrivert.h"
#include "trfelinit.h"
#include "hermtrf.h"

/**
  The function returns required conductivity coefficient
  from the material conductivity %matrix.

  @param phi - relative humidity at given integration point
  @param t - temperature at given integration point
  @param cm - required conductivity %matrix
  @param ri - block row index of transported media
  @param ci - block column index of transported media
  @param im - index of material (default 0 for homogeneous domain)

  @return The function returns required conductivity coefficient.

  Created 04.2011 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
void condmat(double phi, double t, double *cm, long ri, long ci, long im)
{
  matrix d(2,2); // conductivity matrix for 2D problem
  Tm->ip[im].av[0] = phi;
  if (Tp->ntm > 1)
    Tm->ip[im].av[1] = t;
  Tm->matcond(d, im, ri, ci);
  cm[0] = d[0][0];
  cm[1] = d[0][1];
  cm[2] = d[1][0];
  cm[3] = d[1][1];
  return; 
}



/**
  The function returns required capacity coefficient
  from the material capacity %matrix.
  @param phi - relative humidity at given integration point
  @param t - temperature at given integration point
  @param ri - block row index of transported media
  @param ci - block column index of transported media
  @param im - index of material (default 0 for homogeneous domain)

  @return The function returns required capacity coefficient.

  Created 04.2011 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
double capcoeff(double phi, double t, long ri, long ci, long im)
{
  double ret;
  Tm->ip[im].av[0] = phi;
  if (Tp->ntm > 1)
    Tm->ip[im].av[1] = t;
  ret = Tm->capcoeff(im, ri, ci);
  return ret;
}



/** 
  The function initializes the TRFEL from the HERMES code.
  It serves for the separation of codes due to conflicts in the 
  used type names.

  @param argc - number of arguments passed by the command line
                (it has to be 2 or 3)
  @param argv - array of pinters to strings with command line arguments
                argv[0] has to be string with the program name
                argv[1] has to be string with the input file name
                argv[2] is optional and it should contain keyword specificator (see trfelinit.cpp for more details)

  @return The function does not return anything.

  Created 04.2011 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
void hermtrf_init(int argc, const char *argv[])
{
  
  fprintf (stdout,"\n\n ******************************************************\n");
  fprintf (stdout,"            ____ ___ _____ _____ _\n");
  fprintf (stdout,"           / ___|_ _|  ___| ____| |\n");
  fprintf (stdout,"           \\___ \\| || |_  |  _| | |\n");
  fprintf (stdout,"            ___) | ||  _| | |___| |___\n");
  fprintf (stdout,"           |____/___|_|   |_____|_____| HERMTRF");
  fflush(stdout);
  fprintf (stderr,"\n\n ******************************************************\n");

  // initialize global variables to null values
  initnull_globt();
  
  //  program initiation and data reading
  trfel_init (argc,argv);
}



/**
  The function performs a cleanup of the TRFEL global variables
  and it closes opened files.

  Created 04.2011 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
void hermtrf_close()
{
  /*
  delete Tp;
  delete Gtt;
  delete Tt;
  delete Tm;
  delete Tc;
  delete Tb;
  delete Lsrst;
  delete Stt;
  delete Outdt;
  delete Kmat;
  delete Cmat;
  delete Jmat;
  delete Bmat;
  delete Lbt;
  delete Lbat;
  delete Qbt;
  delete Qbat;
  delete Ltt;
  delete Ltat;
  delete Lqt;
  delete Qqt;
  delete Qqat;
  delete Lqat;
  delete Ltett;
  delete Lht;
  delete Qht;
  */
  fclose(Outt);
}
