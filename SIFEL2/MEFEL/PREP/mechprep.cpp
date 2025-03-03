#define EXTERN
#include "prepalias.h"
#include "input.h"
#include "output.h"
#include "globprep.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "mechbclc.h"
#include "outdriverm.h"
#include "globmat.h"
#include "elemhead.h"
#include "globalg.h"
#include "genfile.h"
#include "iotools.h"
#include "kwdset.h"
#include "stacktrace.h"
#include <stdio.h>




/**
  The function checks presence of required sections in the given
  mefel preprocesor file.

  @param in   - pointer to opened input file with problem description

  created by Tomas Koudelka 4.2008, koudelka@cml.fsv.cvut.cz
*/
void check_reqsec(XFILE *in)
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
  descrip d;
  char *logname;
  char *tmp;

  set_prgname(argv[0]);
  fprintf (stdout,"\n\n *** MEFEL PREPROCESSOR ***\n");
  fprintf (stdout," --------------------------\n");

  initnull_glob();
  
  Mp  = new probdesc;
  Gtm = new gtopology;
  Top = new siftop;
  Mt  = new mechtop;
  Mm  = new mechmat;
  Mc  = new mechcrsec;
  Mb  = new mechbclc;
  Outdm = new outdriverm;

  Bar2d = new barel2d;
  Bar3d = new barel3d;
  Beam2d = new beamel2d;
  Beam3d = new beamel3d;
  Spring = new springel;
  Pelt = new planeelemlt;
  Peqt = new planeelemqt;
  Perlt = new planeelemrotlt;
  Pelq = new planeelemlq;
  Peqq = new planeelemqq;
  Perlq = new planeelemrotlq;
  Pesqt = new planeelemsubqt;
  Pqifc = new plquadinterface;
  Cct = new cctelem;
  Dkt = new dktelem;
  Dst = new dstelem;
  Q4pl = new q4plate;
  Asymlq = new axisymlq;
  Asymqq = new axisymqq;
  Asymlt = new axisymlt;
  Shtr = new shelltr;
  Shq = new shellq;
  Ltet = new lintet;
  Lhex = new linhex;
  Qhex = new quadhex;
  Spltr = new soilplatetr;
  Splq = new soilplateq;
  Hexifc = new hexinterface;
  Check_unused_prop = 1;
  
  if (argc < 3){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : mechprep input_file_name output_file name\n\n");
    delete Mp;  delete Mt;  delete Mm;  delete Mc;  delete Mb;
    return(1);
  }
  in = xfopen(argv[1],"r");
  if (in == NULL){
    fprintf (stderr,"\n Input file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    delete Mp;  delete Mt;  delete Mm;  delete Mc;  delete Mb;
    return(2);
  }
  logname = new char[strlen(argv[1])+5];
  strcpy(logname, argv[1]);
  l = long(strlen(logname));
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
  check_reqsec(in);
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
  input_files(in, d);


  // reading of problem description
  xf_setsec(in, bsec_str[begsec_probdesc]);
  fprintf(stdout, "\n\nReading of problem description . . .");
  Mp->read(in);
  
  //  switch of sequential/parallel version
  //  this is the sequential version
  d.paral=0;

  // main input of nodal and element properties
  err = input(in, &d);
  if(err)
  {
    xfclose(in);
    fclose(Log);
    delete Top;    
    return(err);
  }
  
  
  in->kwdmode = sect_mode_seq;
  // reading of output description
  xf_setsec(in, bsec_str[begsec_outdrv]);
  fprintf(stdout, "\n\nReading of outdriver section . . .");
  Outdm->read(in);
  xfclose(in);
  Outdm->conv_sel_prop(Top);
  
  fprintf(stdout, "\n\nOutput of MEFEL input file %s :\n", argv[2]);
  l = long(strlen(argv[2]))+29;
  for (i=0; i<l; i++)
    fprintf(stdout, "-");
  FILE *out = fopen(argv[2], "wt");
  if (out==NULL){
    fprintf (stderr,"\n Output file could not be opened.");
    fprintf (stderr,"\n Try it again!\n\n");
    fclose(Log);
    return(7);
  }

  /*
  // writing of topology
  fprintf(stdout, "\nTopology export . . .");
  Top->exporttop(Gtm);
  fprintf(stdout, " O.K.\n");
  */

  // writing of problem description
  fprintf(stdout, "\nProblem description output . . .");
  Mp->print(out);
  fprintf(stdout, " O.K.\n");

  // main output of nodes, boundary conditions, elements, materials, ...
  output(out, &d);

  // output of time functions
  if (Mp->tprob == growing_mech_structure)
  {  
    fprintf(stdout, "\n\nTime functions output . . .");
    fprintf(out, "# Time functions for dofs and elements\n");
    Gtm->print_gf(out);
    fprintf(stdout, " O.K.\n");
  }
  // output of outdriver data
  fprintf(stdout, "\nOutdriver description output . . .");
  Outdm->print(out);
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

  delete_glob();

  return(0);
}
 
