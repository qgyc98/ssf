#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "galias.h"
#include "iotools.h"
#include "prepalias.h"

char *Fname;
char *Afname;
char *Tmpl;
long Edg;
long Ndom;


void read(XFILE *in)
{
  Fname  = new char[in->give_maxlnsize()+1];
  Afname = new char[in->give_maxlnsize()+1];
  Tmpl   = new char[in->give_maxlnsize()+1];
  memset(Fname, 0, sizeof(*Fname)*(in->give_maxlnsize()+1));
  memset(Fname, 0, sizeof(*Afname)*(in->give_maxlnsize()+1));
  memset(Tmpl, 0, sizeof(*Tmpl)*(in->give_maxlnsize()+1));

  xfscanf(in, " %a", Fname);
  xfscanf(in, " %a", Afname);
  xfscanf(in, " %a", Tmpl);  
  xfscanf(in, "%ld", &Ndom);
  fprintf(stdout, "\nFile name for output : %s\n", Fname);
  fprintf(stdout, "\nNumber of subdomains : %ld\n", Ndom);

  //Edg = 0;
  //fscanf(in, "%ld", &Edg);
}



/**
  The function checks presence of required sections in the given
  mefel preprocesor file.

  @param in   - pointer to opened input file with problem description

  created by Tomas Koudelka 11.2009, koudelka@cml.fsv.cvut.cz
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
}


//JB
int printprep(char *fname, char *tmplfn)
{  
  char *prepfn;
  char *topfn;
  char *icdfn;
  FILE *out;
  XFILE *tmplf;
  long i, k, fmt, ln;
  char str[1205];
  long id_fs, id_pd;
  answertype inicdf_flag;

  tmplf = xfopen(tmplfn,"r");
  if (tmplf==NULL){
    print_err("Template preprocesor file could not be opened.",__FILE__, __LINE__,__func__);
    return(1);
  }
  tmplf->warning = 1;
  tmplf->kwdmode = sect_mode_fwd;
  // detection of sections in the template of preprocesor file
  xfdetect_sect(tmplf, bsec_kwdset, esec_kwdset);
  // checking of sections that have to be detected
  check_reqsec(tmplf);

  ln = strlen(fname);
  prepfn = new char[ln+15];
  topfn  = new char[ln+15];
  icdfn  = new char[ln+15];
  for (i = 0; i < Ndom; i++){
    strcpy(prepfn, fname);
    strcpy(topfn, fname);
    strcpy(icdfn, fname);
    sprintf(prepfn+strlen(prepfn), "%ld.pr", i+1);
    sprintf(topfn+strlen(topfn), "%ld.top", i+1);
    sprintf(icdfn+strlen(icdfn), "%ld.icd", i+1); // initial condition file
    out = fopen(prepfn, "wt");
    if (out == NULL){
      fprintf(stderr, "\nError - unable open preprocessor file %s\n", prepfn);
      fprintf(stderr, " domain number %ld\n", i+1);
      delete [] prepfn;
      delete [] topfn;
      delete [] icdfn;
      return(2);
    }
    // reading of input files 
    xf_setsec(tmplf, bsec_str[begsec_files]);
    id_fs = tmplf->id_sec;
    // topology file name
    fprintf(out, "%s\n", bsec_str[begsec_files].alias);
    fprintf(out, "%s\n", topfn);
    
    // material dbase
    // reading of line with material database file name
    xfscanf(tmplf, " %a", str);
    // write line to the preprocesor file
    fprintf(out, "%s\n", str);
    
    // cross section dbase
    // reading of line with cross-section database file name
    xfscanf(tmplf, " %a", str);
    // write line to the preprocesor file
    fprintf(out, "%s\n", str);
    
    // mesh format and edge/surface numbering
    xfscanf(tmplf, "%k%m", "mesh_format", &meshform_kwdset, &fmt);
    /*
    if (fmt!= t3d)
      {
        sprintf(emsg,"mesh format must be 't3d', see GEFEL/galias.h\n"
		"input file: %s, line: %ld, col: %ld\n", tmplf->fname, tmplf->line, tmplf->col);
	print_err (emsg, __FILE__, __LINE__, __func__);
        return(3);     
      }
*/
    xfscanf(tmplf, "%k%ld", "edge_numbering", &Edg);
    //  mesh format (this is sifel format, see GEFEL/galias.h), edges
    //fprintf(out, "mesh_format t3d\n");
    fprintf(out, "mesh_format %d\n",int(fmt));
    fprintf(out, "edge_numbering %ld\n", Edg);

    // initial condition file
    if (xfscanf(tmplf, "%+k", "inicd_file"))
    {
      fprintf(out, "inicd_file");
      // flag for input file with initial conditions at nodes
      xfscanf(tmplf, "%m", &answertype_kwdset, &inicdf_flag);
      if (inicdf_flag == yes)
      {
        fprintf(out," yes\n%s\n", icdfn);
      }
      else
        fprintf(out," no\n");
    }
    // end of file section  
    fprintf(out, "%s\n\n", esec_str[begsec_files].alias);
    
    // probdesc section
    xf_setsec(tmplf, bsec_str[begsec_probdesc]);
    id_pd = tmplf->id_sec;
    fprintf(out, "%s\n", bsec_str[begsec_probdesc].alias);
    fprintf(out, "domain_number %ld   #domain id", i+1);// domain number
    // do not print \n after comment because th \n must be present after the begsec_probdesc keyword in the template file
    // copy content of the probdesc section
    xf_copysec(tmplf, out);
    // end of probdesc section  
    fprintf(out, "%s\n\n", esec_str[begsec_probdesc].alias);
    
    
    for (k=0; k<tmplf->num_sec; k++)
      {
        if ((k == id_fs) || (k == id_pd))  
	  // these sections have been already copied and modified above
          continue;
        xf_setsec(tmplf, k);
        fprintf(out, "%s",bsec_str[tmplf->asect->name.id].alias);
        // copy section content
        xf_copysec(tmplf, out);
        fprintf(out, "%s\n\n",esec_str[tmplf->asect->name.id].alias);
      }
    fclose(out);
  }
  xfclose(tmplf);
  delete [] prepfn;
  delete [] topfn;
  delete [] icdfn;
  return(0);
}


int main (int argc,char *argv[])
{

  
  //long i, j, ii, jj, k, ie, did;
  XFILE *in;
  
  fprintf (stderr,"\n\n *** PARMEFEL GENERATOR FROM T3D TOPOLOGY ***\n");
  fprintf (stderr," --------------------------------------\n");
  if (argc < 2){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : %s input_file_name\n\n", argv[0]);
    return(1);
  }
  in = xfopen(argv[1],"r");
  if (in==NULL){
    fprintf (stderr,"\n Input file has not been specified or cannot be opened.");
    return(2);
  }
  read(in);
  printprep(Fname, Tmpl); // among others, printprep sets Edg flag used in printtop function
  xfclose(in);
  
  
}
