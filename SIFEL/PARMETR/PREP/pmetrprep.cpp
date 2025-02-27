#define EXTERN

#include <stdio.h>
#include "iotools.h"
#include "globalc.h"
#include "globalg.h"
#include "pglobalc.h"
#include "globprepc.h"
#include "descrip.h"
#include "stacktrace.h"
#include "prepalias.h"



long process_argspc(int argc, const char *argv[]);
void check_reqsec_pc(XFILE *in);
void print_helppc(FILE *out, const char *prgname);




int main(int argc, const char *argv[])
{
  XFILE *in;
  long i, l;
  char *logname, *tmp;
  descrip d;

  Mp    = NULL;
  Tp    = NULL;
  Cp    = NULL;
  Gtu   = NULL;
  Ct    = NULL;
  Cm    = NULL;
  Cmu   = NULL;
  Cml   = NULL;
  Cb    = NULL;
  Outdc = NULL;

  // set name of executable for stacktrace output
  set_prgname(argv[0]);
  fprintf (stdout,"\n\n ******************************************************\n");
  fprintf (stdout,"            ____ ___ _____ _____ _\n");
  fprintf (stdout,"           / ___|_ _|  ___| ____| |\n");
  fprintf (stdout,"           \\___ \\| || |_  |  _| | |\n");
  fprintf (stdout,"            ___) | ||  _| | |___| |___\n");
  fprintf (stdout,"           |____/___|_|   |_____|_____| PMETPREP");
  fflush(stdout);
  fprintf (stderr,"\n\n ******************************************************\n");

  if (process_argspc(argc, argv))
    abort();
  in = xfopen(argv[1],"r");
  if (in == NULL){
    fprintf (stderr,"\n Input file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }
  

  // create inital log file
  logname = new char[strlen(argv[1])+5];
  strcpy(logname, argv[1]);
  l = strlen(logname);
  tmp = strrchr(logname, '.');
  if (tmp)
    sprintf(tmp+1, "plg");    
  else
    sprintf(logname+l, ".plg");
  Log = fopen(logname, "w");

  // dtect sections in the preprocessor file
  in->warning = 1;
  in->kwdmode = sect_mode_seq;
  //in->kwdmode = sequent_mode;
  in->ignorecase = 1;
  // detection of sections in the preprocesor file
  xfdetect_sect(in, bsec_kwdset, esec_kwdset);
  // checking of sections that have to be detected
  check_reqsec_pc(in);

  Mp = new probdesc;
  Tp = new probdesct;
  Cp = new probdescc;
  //Gtu = new gtopology;
  //Top = new siftop;
  Ct = new couptop;
  //Cm = new coupmat;
  //Cmu = new coupmatu;
  //Cml = new coupmatl;
  //Cb = new coupbclc;
  Outdc = new outdriverc;

  // reading of file names of topology, material and cross section
  // xf_setsec(in, bsec_str[begsec_files]);
  
  // reading of problem description
  xf_setsec(in, bsec_str[begsec_probdesc]);
  fprintf(stdout, "\n\nReading of problem description . . .");
  xfscanf (in,"%k%ld","domain_number",&Ndom); // read domain number
  Cp->read(in);

  /*
  // reading of input files 
  xf_setsec(in, bsec_str[begsec_files]);
  // reading of line with topology file name
  xfscanf(in, " %a", d.topf);
  // reading of line with material database file name
  xfscanf(in, " %a", d.matf);
  // reading of line with cross-section database file name
  xfscanf(in, " %a", d.crf);
  // reading of line with topology file format indicator,
  xfscanf(in, "%k%m", "mesh_format", &meshform_kwdset, &d.meshfmt);
  // reading of line with edge number indicator
  xfscanf(in, "%k%ld", "edge_numbering", &d.redgn);
  */

  in->kwdmode = sect_mode_seq;
  // reading of output description
  xf_setsec(in, bsec_str[begsec_outdrv]);
  fprintf(stdout, "\n\nReading of outdriver section . . .");
  Outdc->read(in);
  xfclose(in);

  // main input of nodal and element properties
  //
  // not yet implemented
  /*
  err = inputc(in, &d);
  if(err)
  {
    xfclose(in);
    fclose(Log);
    delete Top;    
    return(err);
  }
  */
  
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

  /*
  // writing of topology
  fprintf(stdout, "\nTopology export . . .");
  Top->exporttop(Gtu);
  fprintf(stdout, " O.K.\n");
  */

  // writing of problem description
  fprintf(stdout, "\nProblem description output . . .");
  fprintf (out,"%d\n",Ndom);
  Cp->print(out);
  fprintf(stdout, " O.K.\n");

  // main output of nodes, boundary conditions, elements, materials, ...
  //
  // to be implemented
  //outputc(out, &d);

  // output of outdriver data
  fprintf(stdout, "\nOutdriver description output . . .");
  Outdc->print(out);
  fprintf(stdout, " O.K.\n");
  fclose(out);

  // close/remove log file
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

  delete Cp;
  delete Ct;
  delete Outdc;

  return 0;
}




/**
  The function processes arguments passed by the command line.

  @param argc[in] - the total number of commandline argumets, i.e. dimension of array argv
  @param argv[in] - pointer to the string array with particular command line arguments

  @ratval 0 - on success
  @retval 1 - error in processing of arguments

  Created by Tomas Koudelka
*/
long process_argspc(int argc, const char *argv[])
{
  switch (argc)
  {
    case 0:
      fflush(stdout);
      fprintf(stderr, "\n\nError - too few arguments are used\n\n");
      print_helppc(stderr, "pmetrprep");
      return 1;
    case 1:
    case 2:
      fflush(stdout);
      fprintf(stderr, "\n\nError - too few arguments are used\n\n");
      print_helppc(stderr,argv[0]);
      return 1;
    case 3:
      break;
    default:
      fflush(stdout);
      fprintf(stderr, "\n\nError - too many arguments is used\n\n");
      print_helppc(stderr,argv[0]);
      return 1;
  }
  return 0;
}



/**
  The function prints help for command line arguments of metr execuatble to the given output text file.
  @param out[in] - pointer to the opened text file
  @param prgname[in] - string with executable name

  @return The function does not return anything. 

  Created by Tomas Koudelka, 08.2021
*/
void print_helppc(FILE *out, const char *prgname)
{
  fflush(out);
  fprintf(out, "\nInput file has not been specified\n");
  fprintf(out, " use : %s infile\n", prgname);
  fprintf(out, " where :\n *   infile is input file name with preprocessor template\n\n");
}



/**
  The function checks presence of required sections in the given
  mefel preprocesor file.

  @param in   - pointer to opened input file with problem description

  created by Tomas Koudelka 8.2021, tomas.koudelka@fsv.cvut.cz
*/
void check_reqsec_pc(XFILE *in)
{
  long err;

  // index of sections must be created
  if (in->index_created != 1)
  {
    print_err("Index of sections has not beeen created", __FILE__, __LINE__, __func__);
    abort();
  }

  /*
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
  */

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
  
  /*
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
  err  = xf_setsec(in, bsec_str[begsec_elvertpr]);
  err += xf_setsec(in, bsec_str[begsec_eledgpr]);
  err += xf_setsec(in, bsec_str[begsec_elsurfpr]);
  err += xf_setsec(in, bsec_str[begsec_elvolpr]);
  if (err == 4)
  {
    print_err("No section with properties of elements has been detected", __FILE__, __LINE__, __func__);
    abort();
  }
  */

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
