#define EXTERN
#include "pglobalt.h"
#include "globalg.h"
#include "prepalias.h"
#include "inputt.h"
#include "outputt.h"
#include "globprept.h"
#include "globalt.h"
#include "globmatt.h"
#include "elemheadt.h"
#include "genfile.h"
#include "iotools.h"
#include "kwdset.h"
#include "stacktrace.h"
#include <stdio.h>

/**
  The function checks presence of required sections in the given
  mefel preprocesor file.

  @param in   - pointer to opened input file with problem description

  created by Tomas Koudelka 9.2010, koudelka@cml.fsv.cvut.cz
*/
void check_reqsec_t(XFILE *in)
{
  long err;

  // index of sections must be created
  if (in->index_created != 1)
  {
    print_err("Index of sections has not beeen created", __FILE__, __LINE__, __func__);
    abort();
  }

  // section with input files must be detected
  err = xf_setsec(in, bsec_str[begsec_files]);  
  switch (err)
  {
    case 1:
      print_err("Section with input files has not been detected", __FILE__, __LINE__, __func__);
      abort();
    default:
      break;
  }

  // section with problem desciption must be detected
  err = xf_setsec(in, bsec_str[begsec_probdesc]);  
  switch (err)
  {
    case 1:
      print_err("Section with problem description has not been detected", __FILE__, __LINE__, __func__);
      abort();
    default:
      break;
  }
  
  // section with description of load cases must be detected
  err = xf_setsec(in, bsec_str[begsec_loadcase]);  
  switch (err)
  {
    case 1:
      print_err("Section with load cases has not been detected", __FILE__, __LINE__, __func__);
      abort();
    default:
      break;
  }

  // at least one section with nodal properties must be detected
  err  = xf_setsec(in, bsec_str[begsec_nodvertpr]);
  err += xf_setsec(in, bsec_str[begsec_nodedgpr]);
  err += xf_setsec(in, bsec_str[begsec_nodsurfpr]);
  err += xf_setsec(in, bsec_str[begsec_nodvolpr]);
  if (err == 4)
  {
    print_err("No section with nodal properties has been detected", __FILE__, __LINE__, __func__);
    abort();
  }

  // at least one section with properties of elements must be detected
  err  = xf_setsec(in, bsec_str[begsec_eledgpr]);
  err += xf_setsec(in, bsec_str[begsec_elsurfpr]);
  err += xf_setsec(in, bsec_str[begsec_elvolpr]);
  if (err == 3)
  {
    print_err("No section with properties of elements has been detected", __FILE__, __LINE__, __func__);
    abort();
  }

  // section with outdriver data must be detected
  err = xf_setsec(in, bsec_str[begsec_outdrv]);
  switch (err)
  {
    case 1:
      print_err("Section with outdriver data has not been detected", __FILE__, __LINE__, __func__);
      abort();
    default:
      break;
  }
}




int main (int argc,char *argv[])
{
  XFILE *in;
  long err, i, l;
  descript d;
  char *logname;
  char *tmp;

  set_prgname(argv[0]);
  fprintf (stdout,"\n\n *** PARTRFEL PREPROCESSOR ***\n");
  fprintf (stdout," -----------------------------\n");


  if (argc < 3){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : ptransprep input_file_name output_file name\n\n");
    return(1);
  }
  in = xfopen(argv[1],"r");
  if (in == NULL){
    fprintf (stderr,"\n Input file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }

  Tp = new probdesct;
  Gtt = new gtopology;
  Top = new siftop;
  Tt = new transtop;
  Tm = new transmat;
  Tc = new transcrsec;
  Tb = new transbclc;
  Outdt = new outdrivert;

  Check_unused_prop = 0;
  
  logname = new char[strlen(argv[1])+5];
  strcpy(logname, argv[1]);
  l = strlen(logname);
  tmp = strrchr(logname, '.');
  if (tmp)
    sprintf(tmp+1, "plg");    
  else
    sprintf(logname+l, ".plg");
  Log = fopen(logname, "w");
  in->warning = 1;
  in->kwdmode = sect_mode_seq;
  //in->kwdmode = sequent_mode;
  in->ignorecase = 1;
  // detection of sections in the preprocesor file
  xfdetect_sect(in, bsec_kwdset, esec_kwdset);
  // checking of sections that have to be detected
  check_reqsec_t(in);
  // checking optional section for materials
  if (xf_setsec(in, bsec_str[begsec_mater]) == 0)
    d.matsec = yes;
  else
    d.matsec = no;

  // checking optional section for cross-section parameters
  if (xf_setsec(in, bsec_str[begsec_crsec]) == 0)
    d.crssec = yes;
  else
    d.crssec = no;

  // reading of input files 
  xf_setsec(in, bsec_str[begsec_files]);
  input_filest(in, d);
  if (d.inicdf == yes)
    Inicdf = d.icf;
  else
    Inicdf = NULL;

  //  creation of object of parallel solver
  Psolt = new psolver (Nproc, Myrank, proc_name, nameLength);
  // reading of problem description
  xf_setsec(in, bsec_str[begsec_probdesc]);
  fprintf(stdout, "\n\nReading of problem description . . .");
  xfscanf (in,"%k%ld","domain_number",&Ndom); // read domain number
  skipline(in); // drop everything until the end of line
  in->kwdmode = sect_mode_fwd;
  if (xfscanf(in,"%+k", "mesprt") == 0) // docasna uprava pro Maderu, mozno odstranit az budou vyresena klicova slova ve funkci Outdt->print
    in->kwdmode = sect_mode_ignore;
  xf_resetsec(in);
  Tp->read(in);
  // reading of basic facts about parallel solver of algebraic equations
  Psolt->read (in,Gtt,Ndom,Tp->tstorkm,Mesprt);

  Lbt   = new linbart;
  Qbt   = new quadbart;
  Ltt   = new trlineart;
  Ltat  = new trlinaxisym;
  Lqt   = new quadlineart;
  Qqt   =  new quadquadrilatt;
  Lqat  = new quadlinaxisym;
  Ltett = new lintett;
  Lht   = new linhext;
  Qht   = new quadhext;
  
  //  switch of sequential/parallel version
  //  this is the parallel version
  d.paral=1;
  if (Psolt->tdd==layered_plate)
    d.paral = 0;
  
  // main input of nodal and element properties
  err = inputt(in, &d);
  if(err)
  {
    xfclose(in);
    fclose(Log);
    delete Top;    
    return(err);
  }

  // reading of general functions
  if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
  {  
    fprintf(stdout, "\n\nReading of section with time functions . . .");
    xf_setsec(in, bsec_str[begsec_loadcase]);
    in->kwdmode = sect_mode_full;
    xfscanf(in, "%k", "time_functions");
    in->kwdmode = sect_mode_seq;
    Gtt->read_gf (in);
  }

  in->kwdmode = sect_mode_seq;
  // reading of output description
  xf_setsec(in, bsec_str[begsec_outdrv]);
  fprintf(stdout, "\n\nReading of outdriver section . . .");
  Outdt->read(in);
  xfclose(in);

  fprintf(stdout, "\n\nOutput of TRFEL input file %s :\n", argv[2]);
  l = strlen(argv[2])+29;
  for (i=0; i<l; i++)
    fprintf(stdout, "-");
  FILE *out = fopen(argv[2], "wt");
  if (out==NULL){
    fprintf (stderr,"\n Output file could not be opened.");
    fprintf (stderr,"\n Try it again!\n\n");
    fclose(Log);
    return(7);
  }

  // writing of topology
  fprintf(stdout, "\nTopology export . . .");
  Top->exporttop(Gtt);
  fprintf(stdout, " O.K.\n");

  // writing of problem description
  fprintf(stdout, "\nProblem description output . . .");
  fprintf (out,"%d\n",Ndom);
  Tp->print(out);
  Psolt->print (out);
  fprintf(stdout, " O.K.\n");

  // main output of nodes, boundary conditions, elements, materials, ...
  outputt(out, &d);

  // output of time functions
  if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
  {  
    fprintf(stdout, "\nTime functions output . . .");
    fprintf(out, "# Time functions for dofs and elements\n");
    Gtt->print_gf(out);
    fprintf(stdout, " O.K.\n");
  }
  // output of outdriver data
  fprintf(stdout, "\nOutdriver description output . . .");
  Outdt->print(out);
  fprintf(stdout, " O.K.\n");
  fclose(out);
  fflush(Log);
  if (ftell(Log) == 0) // no logging was performed
  {
    fclose(Log);
    remove(logname);
  }
  else  // some logging was performed => print message
  {
    fprintf(stdout, "\nSome warnings were registered.\nPlease check the file '%s'\n", logname);
    fclose(Log);
  }
  delete [] logname;
  printf ("\n");
  fprintf(stdout, "\n");
  fprintf(stderr, "\n");



  delete Tt;  
  delete Gtt;  
  delete Top;
  delete Tm;  
  delete Tc;  
  delete Tb;
  delete Outdt;

  delete Lbt;
  delete Qbt;
  delete Ltt;
  delete Ltat;
  delete Lqt;
  delete Qqt;
  delete Lqat;
  delete Ltett;
  delete Lht;
  delete Qht;

  delete Tp;
  delete Psolt;
  return(0);
}
 
