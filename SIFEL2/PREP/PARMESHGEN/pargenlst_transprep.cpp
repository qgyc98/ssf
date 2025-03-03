#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "galias.h"
#include "iotools.h"
#include "prepalias.h"
#include "inputt.h"
#include "descript.h"
#include "globalt.h"
#include "globprept.h"



/**
  The function checks presence of required sections in the given
  mefel preprocesor file.

  @param in   - pointer to opened input file with problem description

  created by Tomas Koudelka 7.2014, koudelka@cml.fsv.cvut.cz
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

  // section with problem description must be detected
  err = xf_setsec(in, bsec_str[begsec_probdesc]);  
  switch (err)
  {
    case 1:
      print_err("Section with problem description has not been detected", __FILE__, __LINE__, __func__);
      abort();
    default:
      break;
  }

  // section with outdriver must be detected
  err = xf_setsec(in, bsec_str[begsec_outdrv]);  
  switch (err)
  {
    case 1:
      print_err("Section with output description has not been detected", __FILE__, __LINE__, __func__);
      abort();
    default:
      break;
  }
}



/**
  The function replaces the last occurence of the asterisk character '*' in the string fname
  with filename by the suffix given in suff. The string have to have allocated sufficient 
  amount of memory to hold additional characters for suffix.

  @param fname - string with filename
  @param suff  - string with required suffix

  @return The function does not return anything but modifies the parameter fname.
 
  Created by Tomas Koudelka, 7.7.2014
*/
void modif_fname(char *fname, char *suff)
{
  char *ar;
  char postf[1000]; 
  long l;

  ar = strrchr(fname, '*');
  postf[0] = '\0';
  postf[999] = '\0';

  if (ar)
  {
    strncpy(postf, ar+1, 999);
    sprintf(ar, "%s", suff);
    l = strlen(suff);
    sprintf(ar+l, "%s", postf);
  }
}



/**
  The function creates according to the preprocesor template file tmplf several preprocessor 
  files with different input/output file names whose suffices are given in the file sufflstf.
  The number of created preprocessor files from the template is given by the number of suffices 
  in the suffix list file. In the template file, each file name containing asterisk character '*'
  is modified so that the asterisk is replaced by the given suffix from the suffix list file.
  
  @param tmplf - pointer to the opened text file with preprocessor template
  @param sufflstf - pointer to the opened text file with the list of file suffices

  @retval 0 - on success
  @retval 1 - some preprocessor file cannot be opened

  Created by Tomas Koudelka, 7.7.2014
*/
int printprep(XFILE *tmplf, XFILE *sufflstf)
{  
  long i, k, nf;
  long id_fs, id_odrv, id_probd;
  char suff[1000];
  char *outfn, *aux, *s;
  descript d;
  FILE *out;

  tmplf->warning = 1;
  tmplf->kwdmode = sect_mode_seq;
  // detection of sections in the template of preprocesor file
  xfdetect_sect(tmplf, bsec_kwdset, esec_kwdset);
  // checking of sections that have to be detected
  check_reqsec(tmplf);
  xf_setsec(tmplf, bsec_str[begsec_files]);
  id_fs = tmplf->id_sec;
  xf_setsec(tmplf, bsec_str[begsec_outdrv]);
  id_odrv = tmplf->id_sec;
  xf_setsec(tmplf, bsec_str[begsec_probdesc]);
  id_probd = tmplf->id_sec;
  tmplf->kwdmode = sect_mode_seq;
  Tp = new probdesct;
  Gtt = new gtopology;
  Tp->read(tmplf);
  Mesprt = 0;
  i = strlen(tmplf->fname);
  outfn = new char[i+1000];
  aux = new char[i+1000];
    
  xfscanf(sufflstf, "%ld", &nf);
  fprintf(stdout, "\n\n Number of processed suffices: %ld\n", nf);
  for (i=0; i<nf; i++)
  {
    // read new suffix
    xfscanf(sufflstf, "%999s", suff);

    // create new preprocesor file name
    aux[0] = '\0';
    strcpy(outfn, tmplf->fname);
    s = strrchr(outfn, '.');
    if (s)
    {
      strcpy(aux,s);
      sprintf(s, "%s", suff);
      sprintf(s+strlen(suff), "%s", aux);
    }
    else
      sprintf(outfn+strlen(outfn), "%s", aux);

    out = fopen(outfn,"wt");
    if (out == NULL){
      print_err("preprocessor file %s cannot be opened.", __FILE__, __LINE__, __func__, outfn);
      return(1);
    }
    else
      fprintf(stdout, " \r  creating %ld. preprocessor file '%s' . . .", i+1, outfn);


    // reading/printing of section with input files 
    xf_setsec(tmplf, bsec_str[begsec_files]);
    input_filest(tmplf, d);
    modif_fname(d.topf, suff);
    if (d.matsec == no)
      modif_fname(d.matf, suff);
    if (d.crssec == no)
      modif_fname(d.crf, suff);
    fprintf(out, "%s\n", bsec_str[begsec_files].alias);
    d.print(out);
    fprintf(out, "%s\n\n", esec_str[begsec_files].alias);

    fprintf(out, "%s\n", bsec_str[begsec_probdesc].alias);
    if (Tp->hdbcont.restore_stat())
      modif_fname(Tp->hdbcont.hdbnamer, suff);
    if (Tp->hdbcont.save_stat())
      modif_fname(Tp->hdbcont.hdbnames, suff);
    Tp->print(out);
    fprintf(out, "%s\n\n", esec_str[begsec_probdesc].alias);
    // reading/printing of remaining sections
    for (k=0; k<tmplf->num_sec; k++)
    {
      if ((k == id_probd) || (k == id_fs) || (k == id_odrv))  
        // these sections have been already/will be copied and modified above/bellow
        continue;
      xf_setsec(tmplf, k);
      fprintf(out, "%s",bsec_str[tmplf->asect->name.id].alias);
      // copy section content
      xf_copysec(tmplf, out);
      fprintf(out, "%s\n\n",esec_str[tmplf->asect->name.id].alias);
    }

    // reading/printing of section with output description
    xf_setsec(tmplf, bsec_str[begsec_outdrv]);
    Outdt->read(tmplf);
    if (Outdt->textout == on)
      modif_fname(Outdt->outfn, suff);
    if (Outdt->gf > grftt_no)
      modif_fname(Outdt->outgrfn, suff);
    if (Outdt->ndiag > 0)
      modif_fname(Outdt->outdiagfn, suff);
    fprintf(out, "%s\n", bsec_str[begsec_outdrv].alias);
    Outdt->print(out);
    fprintf(out, "%s\n\n", esec_str[endsec_outdrv].alias);

    fprintf (stdout, " O.K.");
    fclose(out);
  }
  fprintf(stdout, "\n\n %ld preprocessor files were created sucessfully\n", nf);
  return(0);
}


int main (int argc,char *argv[])
{
  XFILE *tmplf, *sufflstf;
  
  fprintf (stderr,"\n\n *** PARTRFEL GENERATOR FROM TEMPLATE AND DBMAT FILES ***\n");
  fprintf (stderr," --------------------------------------------------------\n");
  if (argc < 3){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : %s prep_template_file_name dbmat_suffix_list_file_name\n\n", argv[0]);
    return(1);
  }

  tmplf = xfopen(argv[1],"r");
  if (tmplf == NULL){
    print_err("template file %s cannot be opened.", __FILE__, __LINE__, __func__, argv[1]);
    return(2);
  }

  sufflstf = xfopen(argv[2],"r");
  if (sufflstf == NULL){
    print_err("dbmat suffices file %s cannot be opened.", __FILE__, __LINE__, __func__, argv[2]);
    return(2);
  }

  Outdt = new outdrivert;
  Tp = new probdesct;
  Mesprt = 0;
  printprep(tmplf, sufflstf);
  
  xfclose(tmplf);
  xfclose(sufflstf);
}
