#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include "iotools.h"
#include "stacktrace.h"



/**
  Function writes error on standard error output device, format specifiers are allowed.

  Parameters:
  @param emsg - string with error message and format specifiers, one terminating character \n will automatically added
  @param errfile - string with source file name where the error was generated
  @param errln - line number in source file  where the error was generated
  @param errfunc - string with function name where the error was generated
  @param ...    - other parameters used in format string emsg passed to the vfprintf function

  @return  Function does not return anything.
  
  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void print_err(const char *emsg, const char *errfile, int errln, const char *errfunc, ...)
{
  va_list params;
  fflush(stdout); // force standard output to be written
  fflush(stderr); // force standard error output to be written
  fprintf(stderr, "\n\nError: ");
  va_start(params, errfunc);
  vfprintf(stderr, emsg, params);
  va_end(params);
  fprintf(stderr, "\n");  
  fprintf(stderr, "detected in source file %s, line %d, function %s\n", errfile, errln, errfunc);
  stack_trace(stderr, 2);
  fflush(stderr);
  fflush(stdout); // force standard output to be written
}



/**
  Function writes warning on standard error output device, format specifiers are allowed in wmsg.

  Parameters:
  @param wmsg - string with warning message, one terminating character \n will automatically added
  @param warnfile - string with source file name where the error was generated
  @param warnln - line number in source file  where the warning was generated
  @param warnfunc - string with function name where the warning was generated
  @param ...    - other parameters used in format string wmsg passed to the vfprintf function

  @return Function does not return anything.
  
  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void print_warning(const char *wmsg, const char *warnfile, int warnln, const char *warnfunc, ...)
{
  va_list params;
  fflush(stdout); // force standard output to be written
  fflush(stderr); // force standard error output to be written
  fprintf(stderr, "\n\nWarning: ");
  va_start(params, warnfunc);
  vfprintf(stderr, wmsg, params);
  va_end(params);
  fprintf(stderr,"\n");
  fprintf(stderr, "detected in file %s, line %d, function %s\n", warnfile, warnln, warnfunc);
  stack_trace(stderr, 2);
  fflush(stderr);
}



/**
  Function writes error on standard error output device.

  Parameters:
  @param n - number of processor
  @param wmsg - string with error message, one terminating character \n will automatically added
  @param errfile - string with source file name where the error was generated
  @param errln - line number in source file  where the error was generated
  @param errfunc - string with function name where the error was generated
  @param ...    - other parameters used in format string emsg passed to the vfprintf function

  @returns Function does not return anything.
  
  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void par_print_err(int n,const char *emsg, const char *errfile, int errln, const char *errfunc, ...)
{
  va_list params;
  fflush(stdout); // force standard output to be written
  fflush(stderr); // force standard error output to be written
  fprintf(stderr, "\n\nError on processor %d:\n", n);
  va_start(params, errfunc);
  vfprintf(stderr, emsg, params);
  va_end(params);
  fprintf(stderr,"\n");
  fprintf(stderr, "detected in file %s, line %d, function %s\n", errfile, errln, errfunc);
  fflush(stderr);
}



/**
  Function writes error on standard error output device.

  Parameters:
  @param n - number of processor
  @param pname - name of processor
  @param emsg - string with error message, one terminating character \n will automatically added
  @param errfile - string with source file name where the error was generated
  @param errln - line number in source file  where the error was generated
  @param errfunc - string with function name where the error was generated
  @param ...    - other parameters used in format string emsg passed to the vfprintf function

  @returns Function does not return anything.
  
  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void par_print_err(int n,const char *pname, const char *emsg, const char *errfile, int errln, const char *errfunc, ...)
{
  va_list params;
  fflush(stdout); // force standard output to be written
  fflush(stderr); // force standard error output to be written
  fprintf(stderr, "\n\nError on processor %d (%s):\n", n, pname);
  va_start(params, errfunc);
  vfprintf(stderr, emsg, params);
  va_end(params);
  fprintf(stderr,"\n");
  fprintf(stderr, "detected in file %s, line %d, function %s\n", errfile, errln, errfunc);
  fflush(stderr);
}



/**
  Function opens file given by parameter name with given mode. This function replaces standard fopen
  function for XFILE type. By default, the kwdmode is set to 'ignore_kwd', warnings are switched off and
  maximum line size is set to 1024.

  Parameters:
  @param name - string with file name
  @param mode - string with file handling mode. Only "r" or "a" reading mode is supported

  @return The function returns pointer to the newly allocated structure XFILE with opened file or NULL
          in case failure of file opening.

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
  modified 11.9.2013, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
XFILE *xfopen(const char *name, const char *mode)
{
  XFILE *ifile;
  if (strchr(mode,'w') || strchr(mode,'a'))
  {
    print_err("unsupported mode for opening of XFILE", __FILE__, __LINE__, __func__);
    abort();
  }
  if (name == NULL)
  {
    return NULL;
  }
  ifile           = new XFILE;
  ifile->file     = fopen(name, mode);
  ifile->line     = 1;
  ifile->col      = 1;
  ifile->lnfpos   = 0;
  ifile->warning  = 0;

  ifile->fname = new char[strlen(name)+1];
  memcpy(ifile->fname, name, strlen(name)+1);

  ifile->kwdmode    = ignore_kwd;
  ifile->ignorecase = 0;
  ifile->set_maxlnsize(1024);
  ifile->maxlnover  = 0;

  ifile->index_created =  0;
  ifile->id_sec        = -1;
  ifile->asect         =  NULL;
  ifile->num_sec       =  0;
  ifile->sect          =  NULL;

  if (ifile->file == NULL)
  {
    print_err("cannot open file %s for reading", __FILE__, __LINE__, __func__, name);
    delete ifile;
    ifile = NULL;
  }
  return ifile;
}



/**
  Function closes file given by parameter f. This function replaces standard fclose
  function for XFILE type.

  Parameters:
  @param f - pointer to opened XFILE structure

  @return The function returns the same value as fclose function

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
  modified 11.9.2013, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
  modified 1.8.2014, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
int xfclose(XFILE* &f)
{
  int ret = 0;
  if (f->maxlnover)
  {
    print_warning("maximum length of line exceeded at least one time\n"
                  "input file: %s, the last oversized line=%ld, col=%ld",
                  __FILE__, __LINE__, __func__, f->fname, f->maxlnover, f->give_maxlnsize());
  }
  if (f->file)
    ret = fclose(f->file);
  f->file = NULL;

  delete f;
  f = NULL;

  return ret;
}



/**
  Function locates the first occurence of keyword given by kwd 
  in the string (line) given by lnstr. Total number of occurences
  is returned via noc.

  Parameters:
  @param lnstr - pointer to scanned string
  @param kwd   - pointer to searched string
  @param ignorecase - flag for case insensitivity (=1)/ case sensitivity (=0)
  @param noc   - number of occurences of searched keyword in lnstr

  @return The function returns either pointer to the first occurence of 
          keyword in the lnstr or NULL in case that the keyword was not found.

  created 02.2008, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
char *locate_kwd(char *lnstr, const char *kwd, long ignorecase, long &noc)
{
  char *ret;
  char *astr;
  size_t kwdl = strlen(kwd);
  size_t strl = strlen(lnstr);
  char *kwd_beg;
  noc = 0;
  ret = NULL;
  astr = lnstr;
  do
  {
    kwd_beg = strstr(astr, kwd, ignorecase);
    if (kwd_beg)
    {
      // testing that the keyword is not part of another word
      if (kwd_beg-lnstr+kwdl < strl)  
      {
        if ((isspace(kwd_beg[kwdl]) == 0) && (ispunct(kwd_beg[kwdl]) == 0))
        // something is joined to the end of located keyword
        {
          astr = kwd_beg+kwdl;
          kwd_beg = NULL;
          continue;
        }
      }
      if (kwd_beg-lnstr > 0)
      {
        if (isspace(kwd_beg[-1]) == 0) 
        // something is prepended to the beginning of located keyword
        {
          if (kwd_beg-lnstr+kwdl < strl)  
          {
            astr = kwd_beg+kwdl;
            kwd_beg = NULL;
            continue;
          }
          else
          {
            kwd_beg = NULL;
            break;
          }
        }
      }
      // located keyword is correct
      if (ret == NULL)
        ret = kwd_beg;
      noc++;
      if (kwd_beg-lnstr+kwdl < strl)  
      // searching continues for other possible keyword occurences
      // in the rest of lnstr
      {
        astr = kwd_beg+kwdl;
        kwd_beg = NULL;
      }
      else
      // end of lnstr was reached -> searching stops
        break;
    }
    else
    // keyword was not found in the lnstr
      break;
  } while(1);

  return ret;
}



/**
  Function detects particular sections in file given by parameter f. Number of 
  searched sections is given by data member n of kwdb and kwde (they must be the same). 
  These structures contain array with base keywords indicating  beginnings or ends 
  of searched sections. Not all of the searched sections have to be detected and 
  the function stores data about detected sections in coresponding data mambers in 
  structure of f (see description of XFILE).

  Parameters:
  @param f - pointer to opened XFILE structure
  @param kwdb    - structure containing array with keywords used for indicating
                   beginning of sections
  @param kwde    - structure containing array with keywords used for indicating
                   end of sections

  @return The function returns number of successfully detected sections

  created  02.2008,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
  modified 11.9.2013, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long xfdetect_sect(XFILE* f, const kwdset &kwdb, const kwdset &kwde)
{
  char *lnstr;
  long i, j;
  size_t l;
  long n = kwdb.n;
  long num_sec;
  long line;
//  long col;
  long lnfpos = 0;
  long lnfpos_p = -1;
  long kwdl;
  long noc;
  char *kwd_beg;
  char *kwd_end;
  xfsection **sect;
  xfsection *tsect;

  switch (f->kwdmode)
  {
    case sect_mode_seq:
    case sect_mode_full:
    case sect_mode_ignore:
    case sect_mode_fwd:
      break;
    default:
      return 0;
  }

  lnstr = f->lnstr;
  sect = new xfsection*[n];
  memset(sect, 0, sizeof(*sect)*n);

  fseek(f->file, 0, SEEK_SET);
  line = 1;
  if (kwdb.n != kwde.n)
  {
    print_err("Number of keywords for beginning of sections is not the same\n"
              " as number of keywords for end of sections\n", __FILE__, __LINE__, __func__);
    abort();
  }
  do
  {
    memset(lnstr, 0, sizeof(*lnstr)*f->give_maxlnsize()+1);
    fscanf(f->file, "%[^\n]",lnstr); //read new line without terminating \n
    clearcomments(lnstr);   // clear possible comments
    // removing trailing \r form the last character of string (due to DOS \r\n newline)
    l = strlen(lnstr);
    if (l > 0)
    {
      if (lnstr[l-1] == '\r')
        lnstr[l-1] = '\0';
    }
    if (strlen(lnstr) == 0) // actual line have zero length
    {
      fscanf(f->file, "%1[\n]", lnstr); // read terminating \n
      line++;
//      col=1;
      lnfpos_p = lnfpos; // backup of previous line beginning position
      lnfpos = ftell(f->file);// backup of beginning line position
      continue;
    }
    for (i=0; i<n; i++)
    {
      // searching beginnings of sections
      kwdl = long(strlen(kwdb.set[i].alias));
      kwd_beg = locate_kwd(lnstr, kwdb.set[i].alias, f->ignorecase, noc);
      if (kwd_beg)
      {
        if ((sect[i] == NULL) && (noc < 2))  // first keyword occurence
        {         
          tsect = new xfsection;
          tsect->beg_secpos = lnfpos + long(kwd_beg - lnstr) + kwdl;
          tsect->beg_secln  = line;
          tsect->beg_seccol = long(kwd_beg - lnstr) + kwdl + 1;
          tsect->name.id = kwdb.set[i].id;
          //tsect->name.alias = new char [kwdl+1];
          //sprintf(tsect->name.alias, "%s",kwdb.set[i].alias);
	  tsect->name.alias = kwdb.set[i].alias;
          tsect->end_secpos = -1;
          sect[i] = tsect;
        }
        else  // multiple keyword occurence
        {
          if (sect[i] == NULL)
            print_err("Multiple occurence of keyword '%s' on line %ld\n", __FILE__, __LINE__, __func__, 
                    kwdb.set[i].alias, line);
          else
            print_err("Multiple occurence of keyword '%s' (line %ld and line %ld)\n", 
                      __FILE__, __LINE__, __func__, kwdb.set[i].alias, sect[i]->beg_secln, line);
          abort();
        }
      }

      // searching ends of sections
      kwdl = long(strlen(kwde.set[i].alias));
      kwd_end = locate_kwd(lnstr, kwde.set[i].alias, f->ignorecase, noc);
      if (kwd_end)
      // keyword for end of section was located
      {
        if (sect[i] == NULL)  
        // error - keyword for end of section preceds keyword for beginning
        {         
          print_err("Keyword for end of section '%s' (line %ld) precedes keyword for beginning of section\n",
                    __FILE__, __LINE__, __func__, kwde.set[i].alias, line);
          abort();
        }
        if ((noc > 1) || (sect[i]->end_secpos >= 0))  // multiple keyword occurence
        {
          if (sect[i]->end_secpos >= 0)
            print_err("Multiple occurence of keyword '%s' on line %ld\n", 
                      __FILE__, __LINE__, __func__, kwde.set[i].alias, line);
          else
            print_err("Multiple occurence of keyword '%s' (line %ld and line %ld)\n", 
                      __FILE__, __LINE__, __func__, kwde.set[i].alias, sect[i]->end_secln, line);

          abort();
        }
        sect[i]->end_secpos = lnfpos + long(kwd_end - lnstr) - 1;
        if ((kwd_end - lnstr) == 0)  // keyword localized at the beginning of new line
        {
          // section ends at the end of previous line
          sect[i]->end_secln  = line-1;
          sect[i]->end_seccol = lnfpos-1-lnfpos_p+1;
        }
        else
        {
          sect[i]->end_secln  = line;
          sect[i]->end_seccol = long(kwd_end - lnstr) + 1;
        }
      }
    }
  } while (!feof(f->file));

  // Counting total number of detected sections
  for(i=0, num_sec=0; i<n; i++)
  {
    if (sect[i])
      num_sec++;
  }
  if (num_sec)
  {
    // Creating index of sections in xfile structure and copying results into it
    f->sect = new xfsection[num_sec];
    j = 0;
    for(i=0; i<n; i++)
    {
      if (sect[i])
      {
        if (sect[i]->end_secpos < 0)
        {
          print_err("Beginning keyword '%s' was found (line %ld) but end keyword '%s' was not found\n",
                    __FILE__, __LINE__, __func__, sect[i]->name.alias, sect[i]->beg_secln, kwde.set[i].alias);
          abort();
        }
        f->sect[j].beg_secpos = sect[i]->beg_secpos; 
        f->sect[j].end_secpos = sect[i]->end_secpos; 
        f->sect[j].beg_secln  = sect[i]->beg_secln; 
        f->sect[j].end_secln  = sect[i]->end_secln; 
        f->sect[j].beg_seccol = sect[i]->beg_seccol; 
        f->sect[j].end_seccol = sect[i]->end_seccol; 
        f->sect[j].name = sect[i]->name;
        j++;
      }
      delete sect[i];
    }
    f->index_created = 1;
    f->id_sec = 0;
    f->asect = &f->sect[0];
    f->num_sec = num_sec;
  }
  delete [] sect;  

  fseek(f->file, f->lnfpos, SEEK_SET);
  return num_sec;
}



/**
  Function resets internal variables of given file so that they 
  points to the beginning of the actual section.

  Parameters:
  @param f - pointer to opened XFILE structure

  @return The function returns nothing

  created  02.2008,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
  modified 11.9.2013, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void xf_resetsec(XFILE *f)
{
  char rfmt[20];   // format string for reading of new line using fscanf
  int fbr;
  long l, errbrk;

  if ((f->index_created == 0) || (f->id_sec < 0) || (f->id_sec > f->num_sec-1))
  {
    if (f->index_created == 0)
      print_err("Error - index of sections in file %s was not created\n", __FILE__, __LINE__, __func__, f->fname);
    if ((f->id_sec < 0) || (f->id_sec > f->num_sec-1))
      print_err("Error - index %ld of section in file %s is out of range (0, %ld)\n", 
                __FILE__, __LINE__, __func__, f->id_sec, f->fname, f->num_sec-1);
     abort();
  }
  fseek(f->file, f->asect->beg_secpos, SEEK_SET);
  f->line   = f->asect->beg_secln;
  f->col    = f->asect->beg_seccol;
  f->lnfpos = f->asect->beg_secpos-f->col+1;

  // fill and prepare line buffer according to xfscanf
  fseek(f->file, f->lnfpos, SEEK_SET);
  sprintf(rfmt, "%%%ld[^\n]%%n", f->give_maxlnsize());
  fscanf(f->file, rfmt, f->lnstr, &fbr); //read new line without terminating \n
  errbrk = check_maxlnsize(f, fbr);
  if (errbrk)
    abort();
  // cut string before end of actual section
  cut_str_sec(f, f->lnstr);
  // removing trailing \r form the last character of string (due to DOS \r\n newline)
  l = long(strlen(f->lnstr));
  if (l > 0)
  {
    if (f->lnstr[l-1] == '\r')
      f->lnstr[l-1] = '\0';
  }
  clearcomments(f->lnstr);   // clear possible comments
  return;
}



/**
  Function sets internal file pointers to the beginning of section with given alias.

  Parameters:
  @param f         - pointer to opened XFILE structure
  @param sec_alias - structure containing string with section name and enum with section alias


  @retval 0 - on success
  @retval 1 - index was not found
  @retval 2 - section index was not built

  created  02.2008,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long xf_setsec(XFILE *f, enumstr sec_alias)
{
  long i, id = -1;

  if (f->index_created == 0)
  // index of sections was not built
    return 2;

  for (i = 0; i < f->num_sec; i++)
  {
    if (strstr(f->sect[i].name.alias, sec_alias.alias) &&
        f->sect[i].name.id == sec_alias.id)
    {
      // required alias of section was found
      id = i;
      break;
    }
  }
  if (id < 0)
  // required alias was not found
    return 1;

  // setup of internal file pointers  
  f->id_sec = id;
  f->asect  = &f->sect[id];
  xf_resetsec(f);
  return 0;
}



/**
  Function sets internal file pointers to the beginning of section of given id

  Parameters:
  @param f      - pointer to opened XFILE structure
  @param sec_id - index of required section in the array sect of XFILE

  @retval 0 - on success
  @retval 1 - index is out of range
  @retval 2 - index was not built

  created  02.2008,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long xf_setsec(XFILE *f, long sec_id)
{
  if (f->index_created == 0)
  // index of sections was not built
    return 2;

  if ((sec_id < 0) || (sec_id > f->num_sec))
  // required index is out of range
    return 1;

  // setup of internal file pointers  
  f->id_sec = sec_id;
  f->asect = &f->sect[sec_id];
  xf_resetsec(f);
  return 0;
}



/**
  Function copises content of the actual section in the file f to the text file out.
  Copying of content starts form the actual file position in the actual section.

  Parameters:
  @param f   - pointer to opened XFILE structure
  @param out - pointer to the opened FILE structure

  @return The function does not return anything

  created  11.2009,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void xf_copysec(XFILE *f, FILE *out)
{
  char *str;
  char fmt[50];
  long pos;
  int size;
  if (f->asect == NULL)
  {
    print_err("cannot determine the actual file section", __FILE__, __LINE__, __func__);
    abort();
  }
  // get current file position
  pos = ftell(f->file);
  if ((pos < f->asect->beg_secpos) || (pos > f->asect->end_secpos))
  {
    // file position is out of the actual section
    print_err("cannot determine position in the actual file section", __FILE__, __LINE__, __func__);
    abort();
  }

  str = new char[1024];
  memset(str, 0, sizeof(*str)*1024);
  // copy body of the section by 1024B blocks
  while (pos+1024 < f->asect->end_secpos)
  {
    pos += 1024;
    fscanf(f->file, "%1024c", str);
    fprintf(out, "%.1024s", str);
  }
  
  // copy rest of the section content
  // size of section content rest in bytes, size is <= 1024
  size = int(f->asect->end_secpos - pos + 1);  
  // create format string "%{n}c" where {n} is replaced by decimal size
  sprintf(fmt, "%%%dc", size);
  fscanf(f->file, fmt, str);  
  fprintf(out, "%.*s", size, str);
  delete [] str;
}



/**
  Function extends standard fscanf function about error checking, input type checking,
  and keywords. Several new conversions are defined or redefined :
  %a - reads string including whitespaces until the newline chracter is reached
       argument of char* type is expected and it has to have sufficient memory allocated to hold
       whole scanned string including trailing \0. 
       *  Format string " %a" skips the initial whitespaces  like 
          ' ', '\t', '\n', '\r' and '\v' and starts the scanning at the first nonwhitespace character
       *  Format string "% a" skips '\n' and '\r' BUT NOT ' ', '\t' and '\v' and than starts scanning. 
          In this case, the scanned string can contain initial spaces or tabelators.
       *  Format string "% #a" skips '\n' and '\r' BUT NOT ' ', '\t' and '\v' and than starts scanning. 
          No more than # characters is read (terminating \0 is not counted) and the scanned string can 
          contain initial spaces or tabs.
  %k - keyword stored in the corresponding argument will be located in the input file. Keyword locating and 
       handling with keyword is controlled by the XFILE kwdmode option. Argument should be string with 
       required keyword, no conversion is performed.Following conversion should perform conversion and 
       assignment of keyword value to the argument which follows argument with keyword string.
  %m - reads enum value given either by string or corresponding integer value. This conversion processes 
       two arguments. The first argument should be pointer to class kwdset, which contains array with 
       enumstr structures containing predefined alias strings and integer values. These aliases 
       (string or integer values) are searched in the input file. Integer value of located alias is assigned 
       to the second argument which is supposed to be a pointer to the integer type.

  Parameters:
  @param in - pointer to the opened XFILE structure
  @param fmt - format string
  @param ... - arguments used for conversions prescribed in the fmt string

  @return returns number of successfully scanned and assigned items
  
  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
  modified 11.9.2013, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
int xfscanf(XFILE *in, const char *fmt, ...)
{
  long  nvar;      // number of required variables by fmt string
  char *types;     // array with identfied conversion types in fmt
  const char **ptypes;   // array with pointers to % before each identfied conversion types in fmt
  long  idt;       // acual id number of processed conversion in types array
  char *subfmt;    // string with whole actually processed conversion
  int br; // total number of read bytes in sscanf call
  int tbr;         // total bytes read in xfscanf call
  va_list params;  // list of passed parameters where the data will be read to
  void *apar;      // pointer to the acutual read parameter
  char *lnstr;     // string with acutal read line
  char *aptr;      // actual position in lnstr where the data are read from
  long  errbrk;    // break flag for error in sscanf
  long  read;      // flag for reading of new line
  long  ret;       // return value from sscanf or get... functions
  long war;        // warning flag
  long readkwdval; // flag for reading keyword value
  kwdset *akwdset; // keyword set for enums 
  long l;          // lenght of string - auxiliary variable
  long start_pos;  // starting position in the file
  char *akwd;      // pointer to string with actual keyword (used in error messages)
  long aline;      // actual line number before reading of keyword (used for error message)
  long acol;       // actual column number before reading of keyword (used for error message)
  char rfmt[20];   // format string for reading of new line using fscanf
  int fbr;         // number of bytes read while reading of new line using fscanf
  int width;       // total number of character read for "%#a" conversion
  long w;          // number of characters read during one pass of %a conversion
  const char *afmt;      // referes to ordinary format string which caused a matching failure (used for proc_ord_fmt return)
  char *aux;       // auxiliary string for format error messages
  char *auxapar;   // auxiliary pointer to the acutual read parameter (used for "%[...]")
  long nskip;      // number of skipped conversions due to optional keywords
  long optkwd;     // indicator of the optional keyword conversion "%+k"
  long optenum;    // indicator of optional enum conversion handling "%+m" (composed integer values are supported)
  long skip_kwdval;// indicator for skipping of the keyword value conversion due to not found of the optional keyword
  long repeat;     // auxiliary flag for skipping of the keyword value due to not found of the optional keyword
  long multiline;  // flag for reading of multiline strings 
  long brmln;      // total number of read bytes in all sscanf call used for multiline strings (%a conversion)

  // detection of total number of required conversions
  nvar = checkfmt(fmt);
  if (nvar == 0)
    return 0;
    
  types  = new char[nvar];
  ptypes = new const char*[nvar];
  subfmt = new char[strlen(fmt)+5];
  //  types  = (char*)alloca(nvar*sizeof(*types));
  //  ptypes = (char**)alloca(nvar*sizeof(*ptypes));
  //  subfmt = (char*)alloca((strlen(fmt)+5)*sizeof(*subfmt));

  lnstr  = in->lnstr;
  aptr   = lnstr + in->col - 1;
  akwd = NULL;
  
  // identification of each conversions and
  // checking of their correct types 
  gettypes(fmt, types, ptypes, nvar);

  va_start(params, fmt);
  idt = 0;
  errbrk = 0;
  readkwdval = 0;
  if ((in->line == 1) && (in->col == 1))
    read = 1; // file is read for the first time
  else
    read = 0;
  width=0;
  nskip = 0;
  skip_kwdval = 0;
  optkwd = 0;
  multiline = 0;
  brmln = 0;
  br = 0;
  tbr = 0;
  start_pos = ftell(in->file);
  // copy format for actual conversion to the subfmt and append %n due to 
  // column counter
  getsubfmt(ptypes[idt], subfmt,1);
  if ((check_asterisk(subfmt)==1) || (types[idt] == ' '))
  // assignment suppression indicator ("%*...") was detected -> no parameter should be passed
    apar = NULL;
  else
  {
    // otherwise selection of newly read parameter
    apar = va_arg(params, void*);
    if (types[idt] == '[') //create zero length string
      *((char*)apar) = '\0';
  }
  auxapar = (char*)apar;

  do
  {
    if (read && readkwdval && (in->kwdmode == line_mode))
    {  
      // reading of new line is required and the keyword search mode is 
      // line_mode => error keyword value was not found in the actual line
      print_err("cannot read value of keyword '%s' (fmt '%s', par id=%ld)\n"
                "input file: %s, line=%ld, col=%ld", __FILE__, __LINE__, __func__, 
                akwd, fmt, idt+1, in->fname, in->line, in->col);
      abort();      
    }
    /*    if (read && readkwdval && ((in->kwdmode == sect_mode_fwd) || (in->kwdmode == sect_mode_seq)) && 
        (in->line == in->asect->end_secln))
    {  
      // reading of new line is required and the keyword search mode is 
      // sect_mode_seq or sect_mode_fwd => error keyword value was not found in the given section
      print_err("cannot read value of keyword '%s' (fmt '%s', par id=%ld)\n"
                "input file: %s, line=%ld, col=%ld", __FILE__, __LINE__, __func__, 
                akwd, fmt, idt+1, in->fname, in->line, in->col);
      abort();      
      }*/
    if (read) // read new line
    {
      memset(lnstr, 0, sizeof(*lnstr)*in->give_maxlnsize()+1);
      checkfeof(in, fmt, idt, ptypes[idt], subfmt, types, apar, __FILE__, __LINE__, __func__);
      checkeos(in, fmt, idt, ptypes[idt], subfmt, types, apar,  __FILE__, __LINE__, __func__);
      // prepare read format string for fscanf function in the form "%x[^\n]", where x is in->maxlnsize
      if (feof(in->file))
      { 
        // This is the case of an ordinary format which begins with space
        // or optional keyword in the "%+k" format.
        // These cases are not catch in the checkfeof/checkeos.
        lnstr[0] = '\0';
      }
      else
      {
        sprintf(rfmt, "%%%ld[^\n]%%n", in->give_maxlnsize()-(in->col-1));
        fbr = 0;
        fscanf(in->file, rfmt, lnstr, &fbr); //read new line without terminating \n
        errbrk = check_maxlnsize(in, fbr);
        if (errbrk)
          break;
        // cut string before end of actual section
        cut_str_sec(in, lnstr);
        // removing trailing \r form the last character of string (due to DOS \r\n newline)
        l = long(strlen(lnstr));
        if (l > 0)
        {
          if (lnstr[l-1] == '\r')
            lnstr[l-1] = '\0';
        }
        clearcomments(lnstr);   // clear possible comments
      }
      aptr = lnstr;       // setup actual position in lnstr
      read = 0;
    }
    if (strlen(aptr) == 0)// actual line have zero length
    {
      if ((types[idt] != 'a') || a_fmt_skips_whitespaces(ptypes[idt]) || multiline) 
      // formats "%a" and "%#a" accept even empty string, the only exception is the "% a" format    
      {
        checkfeof(in, fmt, idt, ptypes[idt], subfmt, types, apar, __FILE__, __LINE__, __func__);
        checkeos(in, fmt, idt, ptypes[idt], subfmt, types, apar, __FILE__, __LINE__, __func__);
        if (types[idt] == '[')
        {
          if (check_newline_fmt(ptypes[idt])) // format "%[...]" accepts even '\n'
          {
            if (apar != NULL) // no suppression of format like "%*[...]"
            {
              *((char*)apar) = '\n';
              apar = (char*)apar + 1;
              *((char*)apar) = '\0';
            }
          }
          else // format "%[...]" does not accept '\n' -> empty string can be matched
          {
            idt++; // index of next conversion
            continue;
          }
        }
        if (feof(in->file))
        {
          // This is the case of an ordinary format which begins with space
          // or optional keyword in the "%+k" format.
          // These cases are not catch in the checkfeof/checkeos.
          // Let the proc_ord_fmt process end of file and let the process getkwd optional keyword
        }
        else
        {
          fscanf(in->file, "%1[\n]", lnstr); // read terminating \n
          in->line++;
          in->col=1;
          in->lnfpos = ftell(in->file);// backup of beginning line position
          read = 1;
          continue;
        }
      }
    }
    // for usual conversions (real and integer), a new line is required in case of empty string
    if (strchr("efgdikmopsux",types[idt]) && isemptystr(aptr, " \t\n\r"))
    {
      checkfeof(in, fmt, idt, ptypes[idt], subfmt, types, apar, __FILE__, __LINE__, __func__);
      checkeos(in, fmt, idt, ptypes[idt], subfmt, types, apar, __FILE__, __LINE__, __func__);
      if (feof(in->file))
      {
        // This is the case of the optional keyword in the "%+k" format which
        // is not catch in the checkfeof/checkeos.
        // Let the process getkwd optional keyword
      }
      else
      {
        fscanf(in->file, "%1[\n]", lnstr);  // read terminating \n
        in->line++;
        in->col=1;
        in->lnfpos = ftell(in->file);// backup of beginning line position
        read = 1;
        continue;
      }
    }
    // for extended string conversion "% a", a new line is required in case of empty string (due to DOS \r\n newline)
    if ((types[idt] == 'a') && (ptypes[idt][1] == ' ') && isemptystr(aptr, "\r"))
    {
      checkfeof(in, fmt, idt, ptypes[idt], subfmt, types, apar, __FILE__, __LINE__, __func__);
      checkeos(in, fmt, idt, ptypes[idt], subfmt, types, apar, __FILE__, __LINE__, __func__);
      fscanf(in->file, "%1[\n]", lnstr);  // read terminating \n
      in->line++;
      in->col=1;
      in->lnfpos = ftell(in->file);// backup of beginning line position
      read = 1;
      continue;
    }
    ret = 0;
    war = 0;
    
    switch (tolower(types[idt]))
    {
      case 'a': // %a conversion = reading of actual line to string
        l = long(strlen(aptr));
        if (multiline == 0)
        {
          ret = sscanf(ptypes[idt]+1, "%d", &width); // scan the width specifier which can be placed 
                                                      // after the '%' character and the possible ' ' modifier
          w = l;
          if ((ret == 1) && (l > width))
            w = width;
        }
        else
        {
          if (width > 0)
          { 
            if (long(br)+l <= width)
              w = l;
            else
            {
              print_err("maximum length of the scanned multiline string exceeded\n"
                        "reading of %%a cannot be continued correctly (fmt '%s', par id=%ld)\n"
                        "input file: %s, line=%ld, col=%ld",
                        __FILE__, __LINE__, __func__, fmt, idt+1, in->fname, in->line, in->col);
              errbrk = 1;
              readkwdval = 0;
              akwd = NULL;
              break;
            }
          }
          else
            w = l;
        }
        if (apar != NULL) // no suppression (%*) in format conversion
        {
          memcpy((char*)(apar)+brmln, aptr, sizeof(char)*w);
          ((char *)apar)[brmln+w] = '\0';
          multiline = detect_multln((char*)(apar));
        }
        else
          multiline = detect_multln(aptr);

        brmln += w;
        br = w;
        if (in->line == in->maxlnover)
        {
          print_err("maximum length of line exceeded\n"
                    "reading of %%a cannot be continued correctly (fmt '%s', par id=%ld)\n"
                    "input file: %s, line=%ld, col=%ld",
                    __FILE__, __LINE__, __func__, fmt, idt+1, in->fname, in->line, in->col);
          errbrk = 1;
        }
        else
          errbrk = 0;
        readkwdval = 0;
        akwd = NULL;
        break;
      case 'k': // %k conversion = locating of keyword depending on the in->kwdmode and in->ignorecase
        // locates keyword from the actual line position, 
        // in the whole line or in the given section
        aline = in->line;
        acol  = in->col;
        ret = getkwd(in, lnstr, aptr, (char *)apar, br); 
        optkwd = (ptypes[idt][1] == '+');       
        if ((ret && (optkwd == 0)) ||    // an error occured and the keyword was compulsory
            ((ret > 1) && (optkwd)))     // other error than "cannot read keyword" and the keyword was optional
        {
          // in case of error reading keyword
          switch (ret)
          {
            case 1:
              print_err("cannot read keyword '%s' (fmt '%s', par id=%ld)\n"
                        "input file: %s, line=%ld, col=%ld", 
                        __FILE__, __LINE__, __func__, (char *)apar, fmt, idt+1, in->fname, aline, acol);
              break;
            case 2:
              print_err("multiple occurence of keyword '%s' on line (fmt '%s', par id=%ld)\n"
                        "input file: %s, line=%ld, col=%ld", 
                        __FILE__, __LINE__, __func__, (char *)apar, fmt, idt+1, in->fname, aline, acol);
              break;
            case 3:
              print_err("zero length of keyword '%s' (fmt '%s', par id=%ld)\n"
                        "input file: %s, line=%ld, col=%ld", 
                        __FILE__, __LINE__, __func__, (char *)apar, fmt, idt+1, in->fname, aline, acol);
              break;
              
          }
          akwd = NULL;
          errbrk = 1;
        }
        else
        {
          if (ret && optkwd) // keyword was not found but the keyword was optional ("%+k" format)
          {
            akwd = NULL;
            skip_kwdval = 1;
            nskip++;  // the conversion %+k was skipped, the value for this conversion is skipped at the 
                      // end of do-while loop
          }
          else
          {
            // setup flag for reading of keyword value
            readkwdval = 1;
            akwd = (char *)apar;
          }
        }
        break;
      case 'm': // %m conversion = reading of enum aliases either by the string alias or by the integer value
        akwdset = (kwdset *)apar; // the first argument is the keyword set
        // selection of newly read enum parameter
        apar = va_arg(params, void*); // the second argument is real 
        optenum = (ptypes[idt][1] == '+');       
        ret = getenum(aptr, akwdset,(int *)apar, br, war, in->kwdmode, in->ignorecase, optenum);
        errbrk = checkenumerr(in, ret, 1, war, akwdset, fmt, subfmt, idt, 
                              ptypes[idt], __FILE__, __LINE__, __func__);
        readkwdval = 0;
        akwd = NULL;
        break;
      case 'f':  // real number conversions %f, %e, %g, %E
      case 'e':
      case 'g':
        war = checkreal(aptr); // checking whether the actual scanned string contains real number
        if (apar == NULL) // no argument will be read (suppressed (%*) or no format conversion)
          ret = sscanf(aptr, subfmt,&br);
        else
        {
          ret = sscanf(aptr, subfmt,apar,&br);
          errbrk = checkscanferr(in, ret, 1, war, "cannot read real number",
                                 fmt, subfmt, idt, ptypes[idt], __FILE__, __LINE__, __func__);
        }
        readkwdval = 0;
        akwd = NULL;
        break;
      case 'd':   // integer number conversions %d, %i, %D
      case 'i':
        war = checkint(aptr); // checking whether the actual scanned string contains integer number
        if (apar == NULL) // no argument will be read (suppressed (%*) or no format conversion)
          ret = sscanf(aptr, subfmt,&br);
        else
        {
          ret = sscanf(aptr, subfmt,apar,&br);
          errbrk = checkscanferr(in, ret, 1, war, "cannot read integer number",
                                 fmt, subfmt, idt, ptypes[idt], __FILE__, __LINE__, __func__);
        }
        readkwdval = 0;
        akwd = NULL;
        break;
      case 'u':  // unsigned integer number conversion %u
        war = checkuint(aptr);  // checking whether the actual scanned string contains unsigned integer number
        if (apar == NULL) // no argument will be read (suppressed (%*) or no format conversion)
          ret = sscanf(aptr, subfmt,&br);
        else
        {
          ret = sscanf(aptr, subfmt,apar,&br);
          errbrk = checkscanferr(in, ret, 1, war, "cannot read unsigned integer number",
                                 fmt, subfmt, idt, ptypes[idt], __FILE__, __LINE__, __func__);
        }
        readkwdval = 0;
        akwd = NULL;
        break;
      case 'o':  // unsigned integer number conversion %o
        war = checkouint(aptr);  // checking whether the actual scanned string contains unsigned integer number in octal notation
        if (apar == NULL) // no argument will be read (suppressed (%*) or no format conversion)
          ret = sscanf(aptr, subfmt,&br);
        else
        {
          ret = sscanf(aptr, subfmt,apar,&br);
          errbrk = checkscanferr(in, ret, 1, war, "cannot read unsigned integer number",
                                 fmt, subfmt, idt, ptypes[idt], __FILE__, __LINE__, __func__);
        }
        readkwdval = 0;
        akwd = NULL;
        break;
      case 'x':  // unsigned integer number conversion %x, %X, %p
      case 'p':
        war = checkxuint(aptr);  // checking whether the actual scanned string contains unsigned integer number in hexadecimal notation
        if (apar == NULL) // no argument will be read (suppressed (%*) or no format conversion)
          ret = sscanf(aptr, subfmt,&br);
        else
        {
          ret = sscanf(aptr, subfmt,apar,&br);
          errbrk = checkscanferr(in, ret, 1, war, "cannot read unsigned integer number",
                                 fmt, subfmt, idt, ptypes[idt], __FILE__, __LINE__, __func__);
        }
        readkwdval = 0;
        akwd = NULL;
        break;
      case 'n':  // number of read characters conversion %n
        *((int *)apar) = tbr;
        errbrk = 0;
        break;
      case ' ':
        if (idt < nvar-1)
          afmt = proc_ord_fmt(in, ptypes[idt], ptypes[idt+1], lnstr, aptr);
        else
          afmt = proc_ord_fmt(in, ptypes[idt], NULL, lnstr, aptr);
        if (afmt)
        {
          aux = new char[strlen(fmt)];
          if (idt<nvar-1)
            strncpy(aux, afmt, ptypes[idt+1]-afmt);
          else
            strcpy(aux, afmt);
          aux[strlen(aux)] = '\0';
          print_err("cannot read required format string '%s'\n"
                    "(format '%s', subformat '%s', parameter id=%ld)\n"
                    "input file: %s, line=%ld, column=%ld", __FILE__, __LINE__, __func__, 
                    aux, fmt, subfmt, idt+1, in->fname, in->line, in->col);
          delete [] aux;
          errbrk = 1;
        }
        else
          errbrk = 0;
        br = 0;    // correct column number and line were managed by process_ord_fmt function
        readkwdval = 0;
        akwd = NULL;
        break;
      default: // default unscpecified conversion which will not be checked
        if (apar == NULL) // no argument will be read (suppressed (%*) or no format conversion)
          ret = sscanf(aptr, subfmt, &br);
        else
        {
          ret = sscanf(aptr, subfmt, apar, &br);
          if (apar == auxapar) // first pass for "%[...]" format
            errbrk = checkscanferr(in, ret, 1, 0, "cannot read required format",
                                   fmt, subfmt, idt, ptypes[idt], __FILE__, __LINE__, __func__);
        }
        readkwdval = 0;
        akwd = NULL;
    }
    if (errbrk) // an error was encountered leave the loop
      break;
    else
    {
      // conversion was successfully scanned and assigned to the argument

      if (((in->kwdmode == sect_mode_seq) || (in->kwdmode == sect_mode_fwd) || (in->kwdmode == sect_mode_full)) && 
          (readkwdval == 1))
        tbr = in->lnfpos - in->asect->beg_secpos + in->col - 1;
      else
      {
        in->col += br; // increase column number with read bytes
        aptr += br;    // change actual position pointer in the scanned line
        tbr = in->lnfpos - start_pos + in->col - 1;
      }

      // catch the case of format "%*[...]" and nothing was read and nothing else 
      // can be read due to feof or eos
      if ((types[idt] == '[') && (apar == NULL) && (auxapar == apar) && (br == 0))
      {
        checkfeof(in, fmt, idt, ptypes[idt], subfmt, types, apar, __FILE__, __LINE__, __func__);
        checkeos(in, fmt, idt, ptypes[idt], subfmt, types, apar,  __FILE__, __LINE__, __func__);
        errbrk = checkscanferr(in, 0, 1, 0, "cannot read required format",
                               fmt, subfmt, idt, ptypes[idt], __FILE__, __LINE__, __func__);
        break;
      }

      if (multiline) // multiline connecting character was detected during the processing of %a format -> 
        continue;    // go to the next line of the input file and perform the conversion again
        

      if ((types[idt] == '[') && (*aptr == '\0') && check_newline_fmt(ptypes[idt]))  // new line chracters are accepted in the format set "%[...]"
      {
        if (apar)      // no assignment suppression (apar != NULL)
          apar = (char*)(apar) + br;  // change the position for further reading
        else
          auxapar += br;
        if (check_feof_eos(in) == 0)  // this condition allows for multiline reading of format "%[...]"
          continue;
      }

      // process remaining ordinary format from the actual conversion
      if (types[idt] != ' ')
      {
        if (idt < nvar-1)
          afmt = proc_ord_fmt(in, ptypes[idt], ptypes[idt+1], lnstr, aptr);
        else
          afmt = proc_ord_fmt(in, ptypes[idt], NULL, lnstr, aptr);
        if (afmt) // matching failure in proc_ord_fmt()
        {         
          aux = new char[strlen(fmt)];
          if (idt<nvar-1)
            strncpy(aux, afmt, ptypes[idt+1]-afmt);
          else 
            strcpy(aux, afmt);
          aux[strlen(aux)] = '\0';
          print_err("cannot read required format string '%s'\n"
                    "(format '%s', subformat '%s', parameter id=%ld)\n"
                    "input file: %s, line=%ld, column=%ld", __FILE__, __LINE__, __func__, 
                    aux, fmt, subfmt, idt+1, in->fname, in->line, in->col);
          delete [] aux;
          errbrk = 1;
          break;
        }
      }

      // reset variables used for reading of  multiline strings
      brmln = 0;
      width = 0;
      multiline = 0;
      
      do
      {
        repeat = 0;
        idt++; // index of next conversion
        if (idt < nvar)
        {
          // copy format for actual conversion to the subfmt and append %n due to 
          // column counter
          getsubfmt(ptypes[idt], subfmt,1);
          if ((check_asterisk(subfmt)==1) || (types[idt] == ' '))
          // assignment suppression indicator ("%*...") was detected -> no parameter should be passed
            apar = NULL;
          else
          {
            // otherwise selection of newly read parameter
            apar = va_arg(params, void*);
            if (types[idt] == '[') //create zero length string
              *((char*)apar) = '\0';
          }
          // auxiliary pointer to the string beginning used for "%[...]" conversion
          auxapar = (char *)apar;

          if (skip_kwdval)
          {
            nskip++;  // skip the value of the not found optional keyword
            skip_kwdval = 0;
            repeat = 1; // next conversion should be processed
          }
        }
      } while (repeat);
      //tbr += br; 
    }
  } while (idt < nvar); // repeat while the index of actual conversion is less then number of required variables in the format string
  va_end(params);
  ret = idt;
  if (nskip)
    ret -= nskip;
    
  if (types[0] == ' ')
    ret--;
  delete [] types;
  delete [] ptypes;
  delete [] subfmt;
  if (errbrk) 
    abort();
    
  return ret;
}



/**
  Function extends standard fprintf function about line prefix, line postfix, error checking, 
  and keywords. ALL ARGUMENTS USED BY CONVERSIONS MUST BE PASSED BY POINTERS.
  Several new conversions are defined or redefined :
  %a - prints string including the whitespaces 
  %[ - prints string including the whitespaces 
  %k - prints keyword stored in the corresponding argument. The keyword printing is controlled by the 
       XFILE kwdmode option. Argument should be string with required keyword. 
  %m - prints enum value given either by string or corresponding integer value. This conversion processes 
       two arguments. The first argument should be pointer to class kwdset, which contains array with 
       enumstr structures containing predefined alias strings and integer values. These aliases 
       (string integer values) are used for the checking of the second parameter and 
       printed to the output file depending on the XFILE kwdmode setting. Integer value of the printed 
       value is contained in the second argument which is supposed to be a pointer to the integer type.

  Parameters:
  @param out - pointer to the opened XFILE structure
  @param fmt - format string
  @param ... - arguments used for conversions prescribed in the fmt string

  @return returns number of successfully scanned and assigned items
  
  created  11.8.2014,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
int xfprintf(XFILE *out, const char *fmt, ...)
{
  long  nvar;        // number of required variables by fmt string
  char *types;       // array with identfied conversion types in fmt
  const char **ptypes;     // array with pointers to % before each identfied conversion types in fmt
  long  idt;         // acual id number of processed conversion in types array
  char *subfmt;      // string with whole actually processed conversion
  char *pconv;       // pointer to the conversion type character in the subfmt
  unsigned int tbw=0;  // total number of written bytes in xfprintf call
  unsigned int bw;   // total number of written bytes in fprintf call
  long optkwd;       // indicator of the optional keyword conversion "%+k"
  kwdset *akwdset;   // keyword set for enums 
  va_list params;    // list of passed parameters where the data will be read to
  void *apar;        // pointer to the acutual printed parameter
  void *apar2;       // pointer to the second printed parameter (used in the %k conversion)
//  int width;         // total number of character read for "%#a" and "%#s" conversions
  char ct;           // lowercase actual conversion type 
  char aux;

  // detection of total number of required conversions
  nvar = checkfmt(fmt);
  if (nvar == 0)
    return 0;

  types  = new char[nvar];
  ptypes = new const char*[nvar];
  subfmt = new char[strlen(fmt)+5];

  // identification of each conversions and
  // checking of their correct types 
  gettypes(fmt, types, ptypes, nvar);

  va_start(params, fmt);
  idt = 0L;
  optkwd = 0L;

  // copy format for actual conversion to the subfmt and append %n due to 
  // column counter
  getsubfmt(ptypes[idt], subfmt, &pconv, 0);
  if (types[idt] == ' ')
  {
    // ordinary format -> no parameter should be passed
    apar = NULL;
  }
  else
  {
    // otherwise selection of newly read parameter
    apar = va_arg(params, void*);
  }
  apar2 = NULL;

  ct = tolower(types[idt]);
  do
  {
    switch (ct)
    {
      case 'a': // %a conversion = print just ordinary string
      case '[': // %[ conversion = print just ordinary string
        if (pconv[0] == '[')
          pconv[1] = '\0';
        pconv[0] = 's';
        bw = fprintf(out->file, subfmt, apar);
        tbw += bw;       // increase counter of the total bytes written
        out->col += bw;
        break;


      case 's': // %s conversion = print string
        aux = *(pconv+1);  // backup character following the conversion type character
        pconv[1] = '\0';   // cut subformat string 
        bw = fprintf(out->file, subfmt, apar);
        tbw += bw;         // increase counter of the total bytes written
        out->col += bw;
        pconv[1] = aux;    // restore character following the conversion type character
        bw = proc_new_ln_fmt(out, pconv+1); // process remaining format string after conversion
        tbw += bw;         // increase counter of the total bytes written
        break;


      case 'k': // %k conversion = locating of keyword depending on the in->kwdmode and in->ignorecase
        apar2 = NULL;
        optkwd = (ptypes[idt][1] == '+');
        if (optkwd)
          apar2 = va_arg(params, void*);
        if (optkwd && (apar2 == NULL)) // optional keyword and no value is given => do not print optional keyword
          break;          
        else
        {
          pconv[0] = 's'; // change %k to %s due to the printing of the given keyword
          if (out->kwdmode > ignore_kwd) // enable keywords in the output
          {
            bw = fprintf(out->file, subfmt, apar);
            tbw += bw;         // increase counter of the total bytes written
            out->col += bw;
          }
        }
        break;


      case 'm': // %m conversion = reading of enum aliases either by the string alias or by the integer value
        akwdset = (kwdset *)apar; // the first argument is the keyword set
        // selection of newly read enum parameter
        apar = va_arg(params, void*); // the second argument is enum 
        if (akwdset->check_int(*((int*)apar))) // check the correct value of the enum
        {
          print_err("integer value %d is not defined in the given kwdset", __FILE__, __LINE__, __func__, *((int*)apar));
          abort();
        }

        if (out->kwdmode > ignore_kwd) // enable keywords in the output
        {
          pconv[0] = 's'; // change %m to %s due to the printing of the given keyword
          bw = fprintf(out->file, subfmt, akwdset->get_str(*((int*)apar)));
          out->col += bw;
          tbw += bw;       // increase counter of the total bytes written
        }
        else  // do not print keywords
        {
          pconv[0] = 'd'; // change %m to %s due to the printing of the given keyword
          bw = fprintf(out->file, subfmt, *((int*)apar));
          out->col += bw;
          tbw += bw;       // increase counter of the total bytes written
        }
        break;


      case 'f':   // real number conversions %f, %e, %g, %E
      case 'e':
      case 'g':
        if (get_modif(subfmt, pconv, "ll", 2))
        {  // long double type conversion
          bw = fprintf(out->file, subfmt, *((long double*)apar));
          tbw += bw;       // increase counter of the total bytes written
          out->col += bw;
          break;
        }
        if (get_modif(subfmt, pconv, "L", 1))
        {  // long double type conversion
          bw = fprintf(out->file, subfmt, *((long double*)apar));
          tbw += bw;       // increase counter of the total bytes written
          out->col += bw;
          break;
        }
        if (get_modif(subfmt, pconv, "l", 1))
        {  // double type conversion
          bw = fprintf(out->file, subfmt, *((double*)apar));
          tbw += bw;       // increase counter of the total bytes written
          out->col += bw;
          break;
        }
        // float type conversion
        bw = fprintf(out->file, subfmt, *((float*)apar));
        tbw += bw;   // increase counter of the total bytes written
        out->col += bw;
        break;


      case 'd':   // integer number conversions %d, %i, %D
      case 'i':
        if (get_modif(subfmt, pconv, "ll", 2))
        {  // long long int type conversion
           // ISO C++ 1998 does not support ‘long long’ => error
          print_err("'ll' modifier required for %d or %i conversion,\n"
                    " but ISO C++ 1998 printf does not support long long type", __FILE__, __LINE__, __func__);
          /*
          bw = fprintf(out->file, subfmt, *((long int*)apar));
          tbw += bw;       // increase counter of the total bytes written
          out->col += bw;*/
          break;
        }
        if (get_modif(subfmt, pconv, "l", 1))
        {  // double type conversion
          bw = fprintf(out->file, subfmt, *((long int*)apar));
          tbw += bw;       // increase counter of the total bytes written
          out->col += bw;
          break;
        }
        if (get_modif(subfmt, pconv, "hh", 2))
        {  // double type conversion
          bw = fprintf(out->file, subfmt, *((signed char*)apar));
          tbw += bw;       // increase counter of the total bytes written
          out->col += bw;
          break;
        }
        if (get_modif(subfmt, pconv, "h", 1))
        {  // long long int type conversion
          bw = fprintf(out->file, subfmt, *((short int*)apar));
          tbw += bw;       // increase counter of the total bytes written
          out->col += bw;
          break;
        }
        // int type conversion
        bw = fprintf(out->file, subfmt, *((int*)apar));
        tbw += bw;   // increase counter of the total bytes written
        out->col += bw;
        break;


      case 'u':   // unsigned integer number conversion %u, %x, %X, %o
      case 'o':
      case 'x':
        if (get_modif(subfmt, pconv, "ll", 2))         
        {  // long long int type conversion
           // ISO C++ 1998 does not support ‘long long’ => error
          print_err("'ll' modifier required for %o, %o or %x conversion,\n"
                    " but ISO C++ 1998 printf does not support long long type", __FILE__, __LINE__, __func__);
          /*
          bw = fprintf(out->file, subfmt, *((unsigned long int*)apar));
          tbw += bw;       // increase counter of the total bytes written
          out->col += bw;*/
          break;
        }
        if (get_modif(subfmt, pconv, "l", 1))
        {  // double type conversion
          bw = fprintf(out->file, subfmt, *((unsigned long int*)apar));
          tbw += bw;       // increase counter of the total bytes written
          out->col += bw;
          break;
        }
        if (get_modif(subfmt, pconv, "hh", 2))
        {  // double type conversion
          bw = fprintf(out->file, subfmt, *((unsigned char*)apar));
          tbw += bw;       // increase counter of the total bytes written
          out->col += bw;
          break;
        }
        if (get_modif(subfmt, pconv, "h", 1))
        {  // long long int type conversion
          bw = fprintf(out->file, subfmt, *((unsigned short int*)apar));
          tbw += bw;       // increase counter of the total bytes written
          out->col += bw;
          break;
        }
        // unsigned int type conversion
        bw = fprintf(out->file, subfmt, *((unsigned int*)apar));
        out->col += bw;
        tbw += bw;   // increase counter of the total bytes written
        break;


      case 'c':
        if (get_modif(subfmt, pconv, "l", 1))
        {  // wide character type conversion
          bw = fprintf(out->file, subfmt, *((wchar_t *)apar));
          tbw += bw;       // increase counter of the total bytes written
          out->col += bw;
          break;
        }
        // character type conversion
        bw = fprintf(out->file, subfmt, *((char *)apar));
        tbw += bw;   // increase counter of the total bytes written
        out->col += bw;
        break;


      case 'p':   // pointer conversion %p
        bw = fprintf(out->file, subfmt, *((void**)apar));
        tbw += bw;       // increase counter of the total bytes written
        out->col += bw;
        break;


      case 'n':  // number of written characters conversion %n
        *((int*)apar) = tbw;
        break;


      case ' ':  // format string without conversion
        bw = proc_new_ln_fmt(out, subfmt);
        tbw += bw;       // increase counter of the total bytes written
        break;


      default: // default unscpecified conversion = error
        print_err("unknown type of conversion is required", __FILE__, __LINE__, __func__);
        abort();
    }

    idt++; // index of next conversion
    
    if (idt < nvar)
    {
      ct = tolower(types[idt]);

      // copy format for actual conversion to the subfmt and append %n due to 
      // column counter
      getsubfmt(ptypes[idt], subfmt, &pconv, 0);

      if (ct == ' ')  // plain format without conversion does not require parameter
        continue;

      // newly printed parameter
      if (optkwd == 0) 
        apar = va_arg(params, void*);
      else // in the case of optional keyword, the next argument has been already obtained as apar2
      {
        if (apar2 == NULL)  // optional keyword value was not defined => go to the next conversion
        {
          idt++;
          apar = va_arg(params, void*);
        }
        else
          apar = apar2;
        optkwd = 0;
      }
    }
  } while (idt < nvar); // repeat while the index of actual conversion is less then number of required variables in the format string

  return tbw;
}



/**
  Function prints given prefix/postfix strings after/before each
  new line character in the fmt string.

  @param out - pointer to the opened XFILE
  @param fmt - pointer to the format processed
  
  @return The function returns the number of written characters to the out file.

  Created by Tomas Koudelka, 11.8.2014, tomas.koudelka@fsv.cvut.cz
*/
long proc_new_ln_fmt(XFILE *out, const char *fmt)
{
  size_t l = strlen(fmt);
  const char *ptrb = fmt;
  const char *ptre = strpbrk(fmt, "\n");
  char *astr = new char[l+1];
  long ret = 0L;
  int bw;

  memset(astr, 0, sizeof(*astr)*(l+1));

  do
  {
    if (ptre == NULL)
    {
      fprintf(out->file, ptrb, NULL);
      ret += long(fmt+l-ptrb);
      out->col += long(fmt+l-ptrb);
      break;
    }
    else
    {
      strncpy(astr, ptrb, ptre-ptrb);
      astr[ptre-ptrb] = '\0';
      fprintf(out->file, astr, NULL);
      ret += long(ptre-ptrb);
      bw = fprintf(out->file, "%s\n", out->lnpostf);
      ret += bw;
      out->line++;
      out->lnfpos = ftell(out->file);
      out->col = 1;
      bw = fprintf(out->file, "%s", out->lnpref);
      out->col += bw;
      ret += bw;
      ptrb = ptre+1;
      ptre = strpbrk(ptrb, "\n");
    }
  } while(1);

  return ret;
}	



/**
  The function checks for the required conversion modifier given by the modif argument
  in the format string fmt.

  @param fmt - format string with only one conversion
  @param pconv - pointer to the conversion character
  @param modif - pointer to a string with the given conversion modifier
  @param modifl - length of the modif string

  @retval 0 - no corresponding modifier was found
  @retval 1 - the given modifier was found

  Created by Tomas Koudelka, 11.8.2014, tomas.koudelka@fsv.cvut.cz
*/
long get_modif(const char *fmt, const char *pconv, const char *modif, int modifl)
{
  const char *aptr;
  aptr = strstr(fmt, modif);
  
  if (aptr) // modifier character was found
  {
    if (aptr < pconv) // modifier was found before conversion character
    {
      if (modifl == 1)
      {
        if (*(aptr-1) == *modif)
          return 0;  // double character modifier was found but single character was modifier required        
      }

      return 1; // single or double modifier character was found according to required modif string
    }
    else
      return 0;   // modifier was found after conversion character => ordinary character in the format string
  }

  return 0; // modifier was not found
}



/**
  The function checks for presence of variable width of string conversion %s in the whole format string.

  @param fmt - a format string
  @param ptypes - array of pointers to tokens of the parsed format string fmt
  @param types  - array of conversion types
  @param nvar   - length of array types and ptypes
  @param fmtl   - length of string fmt
*/
long get_num_args(const char *fmt, char **ptypes, char *types, long nvar, long fmtl)
{
  long narg = 0L; // the total number of expected arguments by the given conversion types
  long i;
//  char *idw = NULL;
  char ct;

  char **vw = new char*[nvar];
  memset(vw, 0, sizeof(*vw)*nvar);

  int *avw = new int[nvar];
  memset(vw, 0, sizeof(*avw)*nvar);

  char *aux;
  char *subfmt = new char[fmtl+1];
  memset(subfmt, 0, sizeof(*subfmt)*(fmtl+1));
  
  for (i=0; i < nvar; i++)
  {
    ct = tolower(types[i]);
    if (ct == ' ')
      continue;
    narg++;
    if (ct == 'm')
      narg++;
    getsubfmt(ptypes[i], subfmt, 0);
    aux = strpbrk(subfmt, "*");
    if (aux)
    {
      vw[i] = aux;
      aux = strpbrk(subfmt, "$");
      if (aux)
      {
        if (sscanf(vw[i]+1, "%d", avw+i) < 1)
        {
          print_err("cannot read index of argument for variable width in the format string:\n"
                    "\"%s\", at position %ld", __FILE__, __LINE__, __func__, fmt, long(vw[i]+1-fmt));
          abort();
        }
      }
      else
        narg++;
    }
  }
  delete [] vw;
  delete [] avw;
  delete [] subfmt;
  return narg;
}



/**
  Function removes commented characters from the string. The '#' character
  is assumed to be beginning of commented part of string.

  Parameters:
  @param lnstr - pointer to opened XFILE structure

  @return The function does not return anything.

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void clearcomments(char *lnstr)
{
  long i = 0;
  while (lnstr[i] != '\0')
  {
    if (lnstr[i] == '#')
    {
      lnstr[i] = '\0';
      return;
    }
    i++;
  }
  return;
}



/**
  Function scans for format conversions in the given string. 
  The number of conversion is increased by one in case that the
  first character of the format string does not contain a conversion
  indicator (%). In such the case the beginning of the format string
  up to the first conversion indicator is considered as an individual conversion
  (It can be used for skipping initial spaces in the input file - adding initial space before
   conversion indicator - " %a", simple keyword management "keyword %ld" etc.; see description 
   of scanf function for handling of nonconversion charcaters.)
  Syntax of format string is the same as in the standard scanf function.

  Parameters:
  @param fmt - pointer to opened XFILE structure

  @return The function returns number of located conversions.

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long checkfmt(const char *fmt)
{
  long i, l, nt, idp;
  l = long(strlen(fmt));
  nt = 0;
  idp = 0;
  for(i=0; i<l; i++)
  {
    if ((fmt[i] == '%') && (idp == 0))
    {
      idp = 1;
      continue;
    }
    if ((fmt[i] == '%') && (idp == 1))
    {
      if (i==1) // initial chracters are considered as an individual conversion
        nt++;
      idp = 0;
      continue;
    }
    if ((i == 0) && (fmt[i] != '%'))  // initial characters before the first % are considered as an individual conversion
      nt++;
    if ((fmt[i] != '%') && (idp == 1))
    {
      nt++;
      idp = 0;
    }
  }
  return nt;
}



/**
  Function assembles one beginning conversion from the format string and copies it to the subfmt string. 
  If the bc flag is set up than the byte count conversion (%n) is appended to subformat string.

  Parameters:
  @param fmt - format string
  @param subfmt - subformat string
  @param bc - byte count flag whether %n is appended to the subformat string (bc=1) or not (bc=0).

  @return The function does not return anything.

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void getsubfmt(const char *fmt, char *subfmt, long bc)
{
  long i, l;
  const char *afmt;
  
  l = long(strlen(fmt));
  i = 0;
  memset(subfmt, 0, sizeof(*subfmt)*(l+5));

  if (*fmt != '%')
  {
    i=1;
    afmt = NULL;
    do
    {
      if (fmt[i] == '%')
      { 
        if (fmt[i+1] == '%')    
          i++;
        else             
          break;
      }
      i++;
    } while(fmt[i]);
  }
  else
   {
    afmt = strpbrk(fmt, "acCdeEfFgGikmnopsSuxX[");  // check for known conversions
    if (afmt != NULL)
    {
      if (*afmt == '[')  // conversion is terminated by ]
      {
        if (afmt[1] == '^') // but if the first character after '^' is the ']' then the']' belongs to the excluded set
          afmt++;
        if (afmt[1] == ']') // but if the first character after '^' or '[' is the ']' then the']' belongs to the excluded/included  set
          afmt++;
        afmt = strchr(afmt+1, ']'); // start searching of ']' after the [ or [] or [^]
        if (afmt == NULL)
        {
          print_err("format %[...] in %s does not contain closing ]", __FILE__, __LINE__, __func__, fmt);
          abort();
        }
      }
      i = long(afmt - fmt) + 1; // number of characters to be copied
    }
    else
    {
      print_err("unknown format specifier is required", __FILE__, __LINE__, __func__);
      abort();
    }
  }
  strncpy(subfmt, fmt, i);
  subfmt[i] = '\0';
  if (bc)
    strcpy(subfmt+i, "%n");
}



/**
  Function assembles one beginning conversion from the format string and copies it to the subfmt string. 
  If the bc flag is set up than the byte count conversion (%n) is appended to subformat string. Additionally,
  it returns pointer to the conversion type character in the subfmt string.

  Parameters:
  @param fmt - format string
  @param subfmt - subformat string (output parameter)
  @param pconv - pointer to the conversion type character in the subfmt string (output parameter)
  @param bc - byte count flag whether %n is appended to the subformat string (bc=1) or not (bc=0).

  @return The function does not return anything but sets arguments subfmt and pconv.

  created  11.8.2014,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void getsubfmt(const char *fmt, char *subfmt, char **pconv, long bc)
{
  long i, l;
  const char *afmt;
  
  l = long(strlen(fmt));
  i = 0;
  memset(subfmt, 0, sizeof(*subfmt)*(l+5));

  if (*fmt != '%')
  {
    i=1;
    afmt = NULL;
    do
    {
      if (fmt[i] == '%')
      { 
        if (fmt[i+1] == '%')    
          i++;
        else             
          break;
      }
      i++;
    } while(fmt[i]);
  }
  else
   {
    afmt = strpbrk(fmt, "acCdeEfFgGikmnopsSuxX[");  // check for known conversions
    if (afmt != NULL)
    {
      if (*afmt == '[')  // conversion is terminated by ]
      {
        if (afmt[1] == '^') // but if the first character after '^' is the ']' then the']' belongs to the excluded set
          afmt++;
        if (afmt[1] == ']') // but if the first character after '^' or '[' is the ']' then the']' belongs to the excluded/included  set
          afmt++;
        afmt = strchr(afmt+1, ']'); // start searching of ']' after the [ or [] or [^]
        if (afmt == NULL)
        {
          print_err("format %[...] in %s does not contain closing ]", __FILE__, __LINE__, __func__, fmt);
          abort();
        }
      }
      i = long(afmt - fmt) + 1; // number of characters to be copied
    }
    else
    {
      print_err("unknown format specifier is required", __FILE__, __LINE__, __func__);
      abort();
    }
  }
  
  strncpy(subfmt, fmt, i);
  subfmt[i] = '\0';
  if (pconv)
  {
    if (*fmt != '%')
      *pconv = subfmt+i-1;
    else
      *pconv = NULL;
  }
  if (bc)
    strcpy(subfmt+i, "%n");
}



/**
  Function extracts particular conversions from the format string fmt and copies their base character to the array types. 
  In addition, the function stores pointers to the beginning of each conversion in the ptypes array.

  Parameters:
  @param fmt    - format string
  @param types  - array with base characters of particular conversions
  @param ptypes - array of pointers to beginning of particular conversions in the format string fmt
  @param nvar   - number of required variables

  @return The function does not return anything.

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void gettypes(const char *fmt, char *types, const char** ptypes, long nvar)
{
  long i, j, l, idp;
  l = long(strlen(fmt)); // length of format string
  idp = 0;         // flag for percent character (%) occurence
  char *subfmt;    // subformat string
  
  subfmt = new char[l+5];
  j=0;  // conversion counter
  for(i=0; i<l; i++)
  {
    if ((fmt[i] != '%') && (i == 0))
    {
      ptypes[j] = fmt;
      j++;
      continue;
    }
    if ((fmt[i] == '%') && (idp == 0)) // check for double percent character - the first occurence of precent character
    {
      idp = 1;
      continue;
    }
    if ((fmt[i] == '%') && (idp == 1)) // check for double percent character - the second occurence of precent character
    {
      if (i==1)
      {
        ptypes[j] = fmt;
        j++;
      }
      idp = 0;
      continue;
    }
    if ((fmt[i] != '%') && (idp == 1)) 
    // check for double percent character - the second occurence of precent character was not found
    // i.e. normal conversion was located
    {
      ptypes[j] = &fmt[i-1]; // pointer to the beginning of conversion (%)
      j++;
      idp = 0;
    }
  }
  for(i=0; i<nvar; i++)
  {
    getsubfmt(ptypes[i], subfmt, 0); // assemble subformat form the format string
    if (*subfmt != '%')
      types[i] = ' ';
    else
    {
      l = long(strlen(subfmt));
      if (subfmt[l-1] == ']')
        types[i] = '[';
      else
        types[i] = subfmt[l-1];
    }
  }
  delete [] subfmt;
}



/**
  Function checks for asterisk chracter which suppresses following conversion assignment

  Parameters:
  @param fmt    - format string

  @retval 0 - no asignment suppression was found
  @retval 1 - assignment suppression was detected

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long check_asterisk(const char *fmt)
{
  long l = long(strlen(fmt));
  long i;
  long fb = 0; // non-whitespace character flag
  for(i=1; i<l; i++)
  {
    if ((fmt[i] == '*') && (fb == 0))
      return 1;
    if ((fmt[i] == '*') && (fb == 1))
    {
      print_err("invalid position of assignment suppression '*' in xfscanf", __FILE__, __LINE__, __func__);
      abort();
    }
    if (fb == 0)
    {
      if (isspace(fmt[i]) != 0)
        fb = 1;
    } 
  }
  return 0;
}



/**
  Function checks maximum length of read line. It is used in xfscanf function
  after new line reading. The function setup flag in the XFILE structure in case 
  the line length was exceeded. 

  Parameters:
  @param in  - pointer to the structure XFILE with opened text file for reading
  @param fbr - number of read characters

  @retval 0 - line size is O.K.
  @retval 1 - line size was exceeded
  @return In addition to that, the function sets flag of parameter 'in' in case the line size was exceeded.

  Created  11.2009,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long check_maxlnsize(XFILE *in, long fbr)
{
  int c;

  // maxlnsize of characters was reached at the actual line
  if (fbr == in->give_maxlnsize())
  {
    c = getc(in->file);   
    if ((c == '\r') || (c == '\n'))
    { 
      // new line character follows -> limit case -> everything is O.K.
      ungetc(c, in->file);  // push back read character
      return 0;
    }
    else    
    { 
      // no new line chracter follows -> line is too long
      ungetc(c, in->file);  // push back read character
      // report warning if the number of previous oversized line differs from the actual one
      if (in->maxlnover != in->line)
      {
        // set oversized line number indicator
        in->maxlnover = in->line;
        if (in->warning)
        {
          if (in->warning == 2) // consider warning as an error
          {
            print_err("maximum length of line exceeded\n"
                      "input file: %s, line=%ld, col=%ld",
                      __FILE__, __LINE__, __func__, in->fname, in->line, in->give_maxlnsize());
            return 1;
          }
          else
            print_warning("maximum length of line exceeded\n"
                          "input file: %s, line=%ld, col=%ld",
                          __FILE__, __LINE__, __func__, in->fname, in->line, in->give_maxlnsize());
        }
      }
    }
  }
  return 0;
}



/**
  Function checks string lnstr for integer value

  @param lnstr - input string

  @retval 0 - only permited characters was found for integer
  @retval 1 - nonpermited character was found 

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long checkint(const char *lnstr)
{
  long i=0, l;
  l = long(strlen(lnstr));
  while (isspace(lnstr[i]) && (i < l)) // search the first non-whitespace character
    i++;
  if (lnstr[i] == '@')
    return 1;
  while((isspace(lnstr[i])==0) && (i < l)) 
  {
    if (strchr("-+0123456789", lnstr[i])) // check each character whether is permitted for integer number
      i++;
    else
    {
      if (lnstr[i] == '@') // if %a conversion follows the '@' may start multiline string =>
        return 0;          // '@' denotes beginning of multiline string => stop checking 
      else
        return 1;
    }
  }
  return 0;
}



/**
  Function checks string lnstr for unsigned integer value

  Parameters:
  @param lnstr - input string

  @retval 0 - only permited characters was found for integer
  @retval 1 - nonpermited character was found 

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long checkuint(const char *lnstr)
{
  long i=0, l;
  l = long(strlen(lnstr));
  while (isspace(lnstr[i]) && (i < l)) // search the first non-whitespace character
    i++;
  if (lnstr[i] == '@')
    return 1;
  while((isspace(lnstr[i])==0) && (i < l))
  {
    if (strchr("+0123456789", lnstr[i])) // check each character whether is permitted for unsigned integer number
      i++;
    else
    {
      if (lnstr[i] == '@') // if %a conversion follows the '@' may start multiline string =>
        return 0;          // '@' denotes beginning of multiline string => stop checking 
      else
        return 1;
    }
  }
  return 0;
}



/**
  Function checks string lnstr for unsigned integer value given by the octal notation.

  Parameters:
  @param lnstr - input string

  @retval 0 - only permited characters was found for integer
  @retval 1 - nonpermited character was found 

  created  8.2014,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long checkouint(const char *lnstr)
{
  long i=0, l;
  l = long(strlen(lnstr));
  while (isspace(lnstr[i]) && (i < l)) // search the first non-whitespace character
    i++;
  if (lnstr[i] == '@')
    return 1;
  while((isspace(lnstr[i])==0) && (i < l))
  {
    if (strchr("+01234567", lnstr[i])) // check each character whether is permitted for unsigned integer number
      i++;
    else
    {
      if (lnstr[i] == '@') // if %a conversion follows the '@' may start multiline string =>
        return 0;          // '@' denotes beginning of multiline string => stop checking 
      else
        return 1;
    }
  }
  return 0;
}



/**
  Function checks string lnstr for unsigned integer value given by the hexadecimal notation.

  Parameters:
  @param lnstr - input string

  @retval 0 - only permited characters was found for integer
  @retval 1 - nonpermited character was found 

  created  8.2014,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long checkxuint(const char *lnstr)
{
  long i=0, l;
  l = long(strlen(lnstr));
  while (isspace(lnstr[i]) && (i < l)) // search the first non-whitespace character
    i++;
  if (lnstr[i] == '@')
    return 1;
  while((isspace(lnstr[i])==0) && (i < l))
  {
    if (strchr("+0123456789abcdefABCDEFxX", lnstr[i])) // check each character whether is permitted for unsigned integer number
      i++;
    else
    {
      if (lnstr[i] == '@') // if %a conversion follows the '@' may start multiline string =>
        return 0;          // '@' denotes beginning of multiline string => stop checking 
      else
        return 1;
    }
  }
  return 0;
}



/**
  Function checks string lnstr for real number value

  Parameters:
  @param lnstr - input string

  @retval 0 - only permited characters was found for integer
  @retval 1 - nonpermited character was found 

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long checkreal(const char *lnstr)
{
  long i=0, l;
  l = long(strlen(lnstr));
  while (isspace(lnstr[i]) && (i < l)) // search the first non-whitespace character
    i++;
  if (lnstr[i] == '@')
    return 1;
  while((isspace(lnstr[i])==0) && (i < l))
  {
    if (strchr("-+.0123456789eE", lnstr[i])) // check each character whether is permitted for real number
      i++;
    else
    {
      if (lnstr[i] == '@') // if %a conversion follows the '@' may start multiline string =>
        return 0;          // '@' denotes beginning of multiline string => stop checking 
      else
        return 1;
    }
  }
  return 0;
}



/**
  Function checks string lnstr contains only set of chracters defined by the string filter
  Usually in this file, the filter contains whitespace characters and the function detects 'empty' string.

  Parameters:
  @param lnstr - input string

  @retval 1 - string contains only characters from the filter set
  @retval 0 - string contains even other characters than the ones from the filter set

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long isemptystr(const char *lnstr,const char *filter)
{
  long i = 0;
  while(lnstr[i] != '\0')
  {
    if (strchr(filter, lnstr[i]) == NULL) // detect whether the actual character is in the filter set
      return 0;
    i++;
  }
  return 1;
}



/**
  Function checks opened file in for end of file (EOF).

  Parameters:
  @param in - pointer to the structure XFILE with opened text file for reading
  Following parameters are used only for error message 
  @param fmt - format string
  @param idt - index of actually processed conversion in fmt
  @param ptype   - pointer to the actual conversion
  @param subfmt  - string for subformat
  @param types   - array of types of particular conversions
  @param apar    - pointer to storage parameter for actual conversion
  @param errfile - string with source file name where the error was detected
  @param errline - line number in source file where the error was detected
  @param errfunc - string with function name where the error was detected

  @return Function does not return anything.

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void checkfeof(XFILE *in, const char *fmt, long idt, const char *ptype, char *subfmt, char *types, void *apar,
               const char *errfile, int errline, const char *errfunc)
{
  // if the the ordinary type is processed and the first character of the format is space
  if ((types[idt] == ' ') && (ptype[0] == ' ')) // then the feof check is performed 
    return;                                       // in the proc_ord_fmt directly (empty string si accepted for the format " ")

  // if the the optional keyword conversion is processed ("%+k")
  if ((types[idt] == 'k') && (ptype[1] == '+')) // let the getkwd skip the keyword and feof will be catch 
    return;                                         // by the next conversion

  if (feof(in->file))
  {    
    getsubfmt(ptype, subfmt, 0);
    if (tolower(types[idt]) == 'k')
      print_err("end of file was reached - reading of keyword '%s' cannot be continued\n"
                "(format '%s', subformat '%s', parameter id=%ld)\n"
                "input file: %s, line=%ld, column=%ld", 
                errfile, errline, errfunc, (char *)apar, fmt, subfmt, idt+1, in->fname, in->line, in->col);
    else 
      print_err("end of file was reached - reading cannot be continued\n"
                "(format '%s', subformat '%s', parameter id=%ld)\n"
                "input file: %s, line=%ld, column=%ld", 
                errfile, errline, errfunc, fmt, subfmt, idt+1, in->fname, in->line, in->col);
    abort();
  }
  return;
}



/**
  Function checks opened file in for end of section 
  (only for sect_mode_seq, sect_mode_fwd, sect_mode_full and sect_mode_ignore).

  Parameters:
  @param in - pointer to the structure XFILE with opened text file for reading
  Following parameters are used only for error message 
  @param fmt - format string
  @param idt - index of actually processed conversion in fmt
  @param ptype   - pointer to the actual conversion
  @param subfmt  - string for subformat
  @param types   - array of types of particular conversions
  @param apar    - pointer to storage parameter for actual conversion
  @param errfile - string with source file name where the error was detected
  @param errline - line number in source file where the error was detected
  @param errfunc - string with function name where the error was detected

  @return Function does not return anything.

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void checkeos(XFILE *in, const char *fmt, long idt, const char *ptype,char *subfmt, char *types, void *apar,
              const char *errfile, int errline, const char *errfunc)
{
  switch (in->kwdmode)
  {
    case ignore_kwd:
    case sequent_mode:
    case line_mode:
      return;
    default:
      break;
  }
  // if the the ordinary type is processed and the first character of the format is space
  if ((types[idt] == ' ') && (ptype[0] == ' ')) // then the eos check is performed 
    return;                                          // in the proc_ord_fmt directly (empty string si accepted for the format " ")

  // if the the optional keyword conversion is processed ("%+k")
  if ((types[idt] == 'k') && (ptype[1] == '+')) // let the getkwd skip the keyword and feof will be catch 
    return;                                         // by the next conversion

  long pos = ftell(in->file);
  if ((in->line == in->asect->end_secln) && 
      (pos > in->asect->end_secpos))
  {
    getsubfmt(ptype, subfmt, 0);
    if ((tolower(types[idt]) == 'k'))
      print_err("end of section '%s' was reached - reading of keyword '%s' cannot be continued\n"
                "(format '%s', subformat '%s', parameter id=%ld)\n"
                "input file: %s, line=%ld, column=%ld",
                errfile, errline, errfunc, in->asect->name.alias, (char *)apar, fmt, subfmt, idt+1, 
                in->fname, in->line, in->col);
    else
      print_err("end of section '%s' was reached - reading cannot be continued\n"
                "(format '%s', subformat '%s', parameter id=%ld)\n"
                "input file: %s, line=%ld, column=%ld",
                errfile, errline, errfunc, in->asect->name.alias, fmt, subfmt, idt+1, 
                in->fname, in->line, in->col);
    abort();
  }
  return;
}



/**
  Function checks opened file in for end of file and end of section 
  (only for sect_mode_seq, sect_mode_fwd, sect_mode_full and sect_mode_ignore).

  Parameters:
  @param in - pointer to the structure XFILE with opened text file for reading

  @retval 0 - no end of file nor end of section were encountered
  @retval 1 - end of file or end of section were encountered

  created  04.2011,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long check_feof_eos(XFILE *in)
{
  if (feof(in->file))
    return 1;

  switch (in->kwdmode)
  {
    case ignore_kwd:
    case sequent_mode:
    case line_mode:
      return 0;
    default:
      break;
  }
  long pos = ftell(in->file);
  if ((in->line == in->asect->end_secln) && 
      (pos > in->asect->end_secpos))
    return 1;

  return 0;
}



/**
  Function cuts string lnstr depending on the end of actually selected section.
  
  Parameters:
  @param in - pointer to the structure XFILE with opened text file for reading
  @param lnstr - input string

  @return Function does not return anything.

  created  2.2008,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void cut_str_sec(XFILE *in, char *lnstr)
{
  switch (in->kwdmode)
  {
    case ignore_kwd:
    case sequent_mode:
    case line_mode:
      return;
    default:
      break;
  }
  if (in->asect->end_secln == in->line)
    lnstr[in->asect->end_seccol-1]='\0';
  return;
}



/**
  Function checks opened file in for error generated by scanf function.

  Parameters:
  @param in - pointer to the structure XFILE with opened text file for reading
  @param ret - returned value from the scanf function
  @param reqret - required returned value from the scanf function which was excepted
  @param war - warning indicator (0=warning ignored/ 1=warning / 2=warning considered as error)
  @param partmsg  - string with part of error message

  Following parameters are used only for error message 
  @param fmt - format string
  @param idt - index of actually processed conversion in fmt
  @param ptype   - pointer to the actual conversion
  @param subfmt  - string for subformat
  @param errfile - string with source file name where the error was detected
  @param errline - line number in source file where the error was detected
  @param errfunc - string with function name where the error was detected

  @retval 0 - no error was encountered
  @retval 1 - error was encountered

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long checkscanferr(XFILE *in, long ret, long reqret, long war, const char *partmsg,
                   const char *fmt, char *subfmt, long idt, const char *ptype, const char *errfile,
                   int errline, const char *errfunc)
{
  if (ret != reqret) // required return value was not obtained
  {
    getsubfmt(ptype, subfmt, 0);
    print_err("%s\n(format '%s', subformat '%s', parameter id=%ld)\n"
              "input file: %s, line=%ld, column=%ld", errfile, errline, errfunc, 
              partmsg, fmt, subfmt, idt+1, in->fname, in->line, in->col);
    return 1;
  }
  if (war && in->warning) // warning was encountered and reporting of warning is switched on
  {
    getsubfmt(ptype, subfmt, 0);
    if (in->warning == 2) // consider warning as an error
    {   
      print_err("%s\n"
                "(format '%s', subformat '%s', parameter id=%ld)\n"
                "input file: %s, line=%ld, column=%ld", errfile, errline, errfunc, 
                partmsg, fmt, subfmt, idt+1, in->fname, in->line, in->col);
      return 1;
    }
    if (in->warning == 1) // only report warning
    {   
      print_warning("%s\n"
                    "(format '%s', subformat '%s', parameter id=%ld)\n"
                    "input file: %s, line=%ld, column=%ld", errfile, errline, errfunc, 
                    partmsg, fmt, subfmt, idt+1, in->fname, in->line, in->col);
    }
  }
  return 0;
}





/**
  Function scans lnstr for keyword. Required keyword is given by apar string. 

  Parameters:
  @param in - pointer to the structure XFILE with opened text file for reading
  @param lnstr - string with line scanned for keyword
  @param aptr  - actual position in lnstr string where should scanning start
  @param apar  - string with required keyword
  @param br    - output parameter in which the read bytes are stored

  @retval 0 - keyword was successfully found
  @retval 1 - keyword was not found
  @retval 2 - multiple keyword detection in lnstr in case line_mode handling
  @retval 3 - invalid keyword is required (zero length of keyword string)

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long getkwd(XFILE *in,char *lnstr, char *&aptr,char *apar, int &br)
{
  long l, tl, lkwd, noc;
  char *kwd_beg;
  char *kwd_beg2;
  char *tmpstr;
  switch (in->kwdmode)
  {
    case ignore_kwd: // ignore keywords
    case sect_mode_ignore:
      br = 0;
      break;
    case sequent_mode: // sequentially scan for keyword
      br = 0;
      lkwd = long(strlen(apar));
      if (lkwd < 1)  // error - zero length of keyword string
        return 3;
      kwd_beg = NULL;
      l = long(strlen(aptr));
      tmpstr = new char[l+1];
      memset(tmpstr, 0, sizeof(*tmpstr)*(l+1));
      sscanf(aptr, "%s", tmpstr);                // selection of the first substring in lnstr
      tl = long(strlen(tmpstr));                        // length of scanned substring, it should be equal to length of keyword string
      if (tl == lkwd)
        kwd_beg = strstr(tmpstr, apar, in->ignorecase); // detection of keyword in the first substring
      delete [] tmpstr;
      if (kwd_beg == NULL)
        return 1;
      else
      {
        kwd_beg = strstr(aptr, apar, in->ignorecase); // detection of keyword in the lnstr -> br
        br = long(kwd_beg-aptr)+lkwd;
      }
      break;
    case line_mode:
      br = 0;
      lkwd = long(strlen(apar));
      if (lkwd < 1)  // error - zero length of keyword string
        return 3;
      kwd_beg = strstr(lnstr, apar, in->ignorecase); // detection of keyword
      if (kwd_beg == NULL)
        return 1;
      l = long(strlen(lnstr));
      tmpstr = new char[l+1];
      memset(tmpstr, 0, sizeof(*tmpstr)*(l+1));
      sscanf(kwd_beg, "%s", tmpstr);                // selection of the first substring in lnstr
      tl = long(strlen(tmpstr));                        // length of scanned substring, it should be equal to length of keyword string
      if (tl != lkwd)
      {
        delete [] tmpstr;
        return 1;
      }
      // detection of multiple occurence of keyword
      kwd_beg2 = strstr(kwd_beg+lkwd, apar, in->ignorecase);
      if (kwd_beg2) 
      {
        memset(tmpstr, 0, sizeof(*tmpstr)*(l+1));
        sscanf(kwd_beg2, "%s", tmpstr);             // selection of the first substring in lnstr
        tl = long(strlen(tmpstr));                        // length of scanned substring, it should be equal to length of keyword string
        if (tl != lkwd)
        {
          delete [] tmpstr;
          return 2;
        }
      }
      delete [] tmpstr;
      br = int(kwd_beg - lnstr) + int(lkwd);
      break;
    case sect_mode_seq:
    case sect_mode_fwd:
    case sect_mode_full:
      br = 0;
      lkwd = long(strlen(apar));
      if (lkwd < 1)  // error - zero length of keyword string
        return 3;
      noc = getkwd_sect(in, lnstr, aptr, apar, 0);
      if (noc == 0) 
        return 1;
      if (noc > 1)
        return 2;
      break;      
    default:
      print_err("Unknown type of keyword handling mode is required\n", __FILE__, __LINE__, __func__);
      abort();
  }
  return 0;
}



/**
  Function scans acutal section of file in for given keyword
  and return number of located keywords. If the lnstr == NULL,
  the function is supposed tobe called outside from the xfscanf,
  otherwise it is supposed to be called by the xfscanf. 
  At the function entry it is supposed:
  - in->line in->col have to be correct actual line and column number
  - in->lnfpos has to be correct file position of the actual line beginning
  - file poition (ftell) corresponds to the last character on the actual line before \n
  if called from xfscanf (lnstr != NULL):
  - content of lnstr may be changed so that it contains line with the first keyword occurence
  - aptr is changed so that it points OVER the first located keyword
  - in->line, in->col are changed so that it points OVER the first located keyword
  - in->lnfpos file position of the actual line beginning
  - position in file is set (fseek) to the end of actual line before \n
  if called externally (lnstr == NULL):
  - lnstr and aptr are unchanged
  - in->line, in->col are changed so that it points TO THE FIRST character of first located keyword
  - in->lnfpos file position of the actual line beginning
  - position in file is set (fseek) to the end of actual line before \n


  Parameters:
  @param in - pointer to the structure XFILE with opened text file for reading
  @param lnstr - string scanned for keyword used in the xfscanf
  @param aptr  - actual position in lnstr string where scanning should start
  @param apar  - string with required keyword
  @param detect_allkwd - flag for detection of all keyword occurences(=1) or 
                         only the firstkeyword occurence (=0).
 
  @return number of suceffully located keywords

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
  modified 11.9. 2013,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz

*/
long getkwd_sect(XFILE *in, char *lnstr, char *&aptr, const char * apar, long detect_allkwd)
{
  switch (in->kwdmode)
  {
    case sect_mode_full:
      // for this mode, searching starts always from the beginning of the section
      xf_resetsec(in); 
      break;
    case sect_mode_seq:
    case sect_mode_fwd:
      break;
    default:
      print_err("Unknown type of keyword handling mode is required\n", __FILE__, __LINE__, __func__);
      abort();
  }

  if (in->index_created == 0)
  {
    print_err("Section keyword mode is required but section index is not built\n"
              "Use xfdetect_sect function before calling xfscanf function\n", __FILE__, __LINE__, __func__);
    abort();
  }
  if ((in->id_sec < 0) || (in->id_sec > in->num_sec))
  {
    print_err("Section id is out of range\n", __FILE__, __LINE__, __func__);
    abort();
  }  

  long tnoc = 0;   // total number of keyword occurences
  long noc  = 0;   // number of keyword occurences in one string
  char *kwd_loc = NULL;       // pointer to the first keyword occurence
  size_t kwdl = strlen(apar); // length of keyword
  long l;           // length of astr(line)
  char *loc = NULL; // pointer to the other keyword occurences
  char *astr;       // string with line which is read actually
  char rfmt[20];    // format string for reading of new line using fscanf 
  long br = 0;      // flag for read bytes;
  long fkwd_line;   // line number with the first keyword occurence
  long fkwd_col;    // column number with the first keyword occurence
  long fkwd_lnfpos; // position of the fkwd_line beginning in the file 
  long fkwd_endlnfpos; // position of the fkwd_line end in the file 

  astr = new char[in->give_maxlnsize()+1];
  astr[0] = '\0';
  if ((lnstr == NULL) || (in->kwdmode == sect_mode_full))
  {
    // getkwd_sect is called out of xfscanf function or 
    // for sect_mode_full, the first line of the section has been already read
    strcpy(astr, in->lnstr+in->col-1);
  }
  else  // getkwd_sect is called from the xfscanf function
  {
    // astr should begin on the actual position on the given line
    strcpy(astr, aptr);
  }
  do 
  {
    if (kwd_loc == NULL) // no keyword was located
      kwd_loc = locate_kwd(astr, apar, in->ignorecase, noc);
    else // at least one occurence of keyword was located and multiple occerences are required to be detected
      loc = locate_kwd(astr, apar, in->ignorecase, noc);

    if (kwd_loc && (tnoc==0)) // the first keyword occurence is detected
    {
      if (br && lnstr)
        strcpy(lnstr, astr);
      if (lnstr == NULL)  // fill the line buffer for external call of the function
        strcpy(in->lnstr, astr);

      fkwd_col    = in->col + long(kwd_loc - astr);
      fkwd_line   = in->line;
      fkwd_lnfpos = in->lnfpos;
      fkwd_endlnfpos = ftell(in->file);;
      if (lnstr)
      {
        // for internal call from the xfscanf the position over the first 
        // located keyword is required
        fkwd_col += long(kwdl);
        if (br == 0)
          aptr += kwd_loc - astr + kwdl;
        else
          aptr = lnstr + fkwd_col - 1;
      }
      if (detect_allkwd == 0)
      {
        tnoc = 1;  // one line (astr) may contain multiple keywords 
                   // (they are all detected by locate_kwd) but the first keyword is required only
        break;
      }
      else
        tnoc += noc;
    }
    else  // keyword was not detected or other occurences of keyword are searched
    {
      // sequential keyword searching required and the keyword was not found
      if ((in->kwdmode == sect_mode_seq) && (kwd_loc == NULL))
      {
        // the actual line is non-empty => stop searching and return because keyword was not found
        if (isemptystr(astr, " \t\n\r") == 0)
          break;
      }
      if (loc)  // other occurence of keyword was detected
        tnoc += noc;
    }
    loc = NULL;
    // reading of new line
    if (in->line < in->asect->end_secln)
    {          
      astr[0] = '\0'; // setting string to zero length
      fscanf(in->file, "%*1[\n]"); //read new line terminating character
      in->line++;
      in->col = 1;
      in->lnfpos = ftell(in->file);
      // prepare read format string for fscanf function in the form "%x[^\n]", where x is in->maxlnsize
      sprintf(rfmt, "%%%ld[^\n]", in->give_maxlnsize());
      fscanf(in->file, rfmt, astr); //read new line
      br = 1;
      // cut string before end of actual section
      cut_str_sec(in, astr);
      // removing trailing \r form the last character of string (due to DOS \r\n newline)
      l = long(strlen(astr));
      if (l > 0)
      {
        if (astr[l-1] == '\r')
          astr[l-1] = '\0';
      }
      clearcomments(astr);   // clear possible comments
    }
    else
      break;

  } while (1); 

  if (kwd_loc)  // setting of file position to the end of line with the first located keyword
  {
    in->line   = fkwd_line;
    in->col    = fkwd_col;
    in->lnfpos = fkwd_lnfpos;
    fseek(in->file, fkwd_endlnfpos, SEEK_SET);
  }
  delete [] astr;
  return tnoc;
}




/**
  Function scans lnstr for keyword. Required keyword is given by apar string. 

  Parameters:
  @param in - pointer to the structure XFILE with opened text file for reading
  @param lnstr - string with line scanned for keyword
  @param aptr  - actual position in lnstr string where should scanning start
  @param apar  - string with required keyword
  @param br    - output parameter in which the read bytes are stored

  @retval 0 - keyword was successfully found
  @retval 1 - keyword was not found
  @retval 2 - multiple keyword detection in lnstr in case line_mode handling
  @retval 3 - invalid keyword is required (zero length of keyword string)

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long getkwd_opt(XFILE *in,char *lnstr, char *&aptr,char *apar, int &br)
{
  long l, tl, lkwd, noc;
  char *kwd_beg;
  char *kwd_beg2;
  char *tmpstr;
  switch (in->kwdmode)
  {
    case ignore_kwd: // ignore keywords
    case sect_mode_ignore:
      br = 0;
      break;
    case sequent_mode: // sequentially scan for keyword
      lkwd = long(strlen(apar));
      if (lkwd < 1)  // error - zero length of keyword string
        return 3;
      kwd_beg = NULL;
      l = long(strlen(aptr));
      tmpstr = new char[l+1];
      memset(tmpstr, 0, sizeof(*tmpstr)*(l+1));
      sscanf(aptr, "%s", tmpstr);                // selection of the first substring in lnstr
      tl = long(strlen(tmpstr));                        // length of scanned substring, it should be equal to length of keyword string
      if (tl == lkwd)
        kwd_beg = strstr(tmpstr, apar, in->ignorecase); // detection of keyword in the first substring
      delete [] tmpstr;
      if (kwd_beg == NULL)
        return 1;
      else
      {
        kwd_beg = strstr(aptr, apar, in->ignorecase); // detection of keyword in the lnstr -> br
        br = int(kwd_beg-aptr)+int(lkwd);
      }
      break;
    case line_mode:
      lkwd = long(strlen(apar));
      if (lkwd < 1)  // error - zero length of keyword string
        return 3;
      kwd_beg = strstr(lnstr, apar, in->ignorecase); // detection of keyword
      if (kwd_beg == NULL)
        return 1;
      l = long(strlen(lnstr));
      tmpstr = new char[l+1];
      memset(tmpstr, 0, sizeof(*tmpstr)*(l+1));
      sscanf(kwd_beg, "%s", tmpstr);                // selection of the first substring in lnstr
      tl = long(strlen(tmpstr));                        // length of scanned substring, it should be equal to length of keyword string
      if (tl != lkwd)
      {
        delete [] tmpstr;
        return 1;
      }
      // detection of multiple occurence of keyword
      kwd_beg2 = strstr(kwd_beg+lkwd, apar, in->ignorecase);
      if (kwd_beg2) 
      {
        memset(tmpstr, 0, sizeof(*tmpstr)*(l+1));
        sscanf(kwd_beg2, "%s", tmpstr);             // selection of the first substring in lnstr
        tl = long(strlen(tmpstr));                        // length of scanned substring, it should be equal to length of keyword string
        if (tl != lkwd)
        {
          delete [] tmpstr;
          return 2;
        }
      }
      delete [] tmpstr;
      br = int(kwd_beg - lnstr) + int(lkwd);
      break;
    case sect_mode_seq:
    case sect_mode_fwd:
    case sect_mode_full:
      lkwd = long(strlen(apar));
      if (lkwd < 1)  // error - zero length of keyword string
        return 3;
      noc = getkwd_sect(in, lnstr, aptr, apar, 0);
      if (noc == 0) 
        return 1;
      if (noc > 1)
        return 2;
      break;      
    default:
      print_err("Unknown type of keyword handling mode is required\n", __FILE__, __LINE__, __func__);
      abort();
  }
  return 0;
}



/**
  Function scans lnstr for enum type. Possible aliases for given enum are stored in the akwdset. 

  Parameters:
  @param lnstr - string scanned for keyword
  @param akwdset - pointer to class containing keyword set and set of corresponding integer values which are perimtted for given enum
  @param apar  - pointer to the integer argument to which the detected enum integer value will be assigned
  @param br    - output parameter which the read bytes are stored in
  @param war   - output parameter which the warning flag is stored in
  @param handling - this parameter controls scanning for keyword
                    = 0 - ignore_kwd akwdset -> do not scan for aliases keywords
                    != 0 - scan either for possible keywords or for the integer id
  @param ignorecase - flag for case insensitivity (=1)/ case sensitivity (=0)
  @param optenum - this parameter controls optional handling with enum values
                   = 0 (default) - enum value is scanned and tested whether it match the enum definition in akwdset 
                   != 0 - enum value is considered to be rather composed bit array obtained as a result of
                          biwise addition of particular values from the enum definition
                          In this case, single values corresponding to enum definition can be given either by keyword or 
                          by integer value but composed values can be given by integers (results by bitwise addition) only.
  @retval 1 - enum was successfully scanned
  @retval 2 - cannot read integer value of enum
  @retval 3 - scanned integer value of enum was not found in akwdset
  @retval 4 - nor enum keyword nor integer value of enum was read

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long getenum(const char *lnstr, kwdset *akwdset, int *apar, int &br, long &war, int handling, long ignorecase, long optenum)
{
  long l, n, tbr;
  char *ekwd;
  long ret = 1;
  int  aux;

  war = 0;
  if (apar == NULL) // no suppression (%*) in format conversion
    apar = &aux;
  switch (handling)
  {
    case ignore_kwd:
      // try to read int value of enum 
      war = checkint(lnstr);
      ret = sscanf(lnstr, "%d%n",apar,&br);
      if (ret != 1)
        ret = 2; // cannot read integer value of enum
      else
      {
        if (akwdset)
        {
          if (optenum)
          { 
            if (matchoptekwdint(apar, akwdset))
              ret = 1; // O.K. composed integer value was found and checked
            else
              ret = 3; // composed enum integer value was not found in kwdset
          }
          else
          {
            if (matchekwdint(apar, akwdset))
              ret = 1; // O.K. integer value was found and checked
            else
              ret = 3; // enum integer value was not found in kwdset
          }
        }
      }
      break;
    default:
      // For all cases of keyword handling, the enum keyword have to be read by the sequential mode
      if (akwdset == NULL)
        return 5; // the keyword set was not passed
      l = long(strlen(lnstr));
      ekwd = new char[l+1];
      n = sscanf(lnstr, "%s%n", ekwd, &br);
      l = long(strlen(ekwd));
      tbr = br-l; // number of whitespaces suppressed in above sscanf %s
      if (n != 1)
      {
        ret=4; // cannot read enum
      }
      else
      {
        if (matchekwd(ekwd, akwdset, apar, ignorecase))
          ret = 1; // O.K. kwd was found and checked
        else
        {
          // try to read int value of enum 
          war = checkint(ekwd);
          ret = sscanf(ekwd, "%d%n",apar,&br);
          br += tbr;
          if (ret != 1)
            ret = 4; // nor enum keyword nor integer value of enum was read
          else
          {
            if (optenum)
            { 
              if (matchoptekwdint(apar, akwdset))
                ret = 1; // O.K. composed integer value was found and checked
              else
                ret = 3; // composed enum integer value was not found in kwdset
            }
            else
            {
              if (matchekwdint(apar, akwdset))
                ret = 1; // O.K. integer value was found and checked
              else
                ret = 3; // enum integer value was not found in kwdset
            }
          }
        }
      }
      delete [] ekwd;
  }
  return ret;
}




/**
  Function searches aliases keywords in the string. (used for enum conversion)
  
  Parameters:
  @param ekwd - string with required keyword 
  @param akwdset - pointer to class which contains set of aliases keywords and corresponding integer values
  @param apar - output parameter (pointer to integer) to which the integer value of detected enum keyword will be assigned
  @param ignorecase - flag for case insensitivity (=1)/ case sensitivity (=0)

  @retval 0 - enum keyword was not found
  @retval 1 - enum keyword was successfully scanned

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long matchekwd(const char *ekwd, kwdset *akwdset, int *apar, long ignorecase)
{
  long i, lkwd, tl;
  const char *kwd_beg; 
  long ret=0;
  long pos_ok;
 
  for(i=0; i<akwdset->n; i++)
  {
    kwd_beg = strstr(ekwd, akwdset->set[i].alias, ignorecase);
    if (kwd_beg)
    {
      tl   = long(strlen(kwd_beg));
      lkwd = long(strlen(akwdset->set[i].alias));
      
      // check position of found keyword string 
      // (corectly found keyword cannot be part of another nonwhitespace part of searched string)
      pos_ok = 0;
      if (kwd_beg != ekwd)
      {
        if (isspace(*(kwd_beg-1))) // only space can precede of the found keyword 
          pos_ok = 1;
      }
      else
        pos_ok = 1;

      if ((tl == lkwd) && pos_ok)
      {
        *apar = akwdset->set[i].id;
        ret   = 1;
        break;
      }
    }
  }
  return ret;
}



/**
  Function searches for the given integer value apar in alias integers. (used for enum conversion)
  
  Parameters:
  @param apar - searched integer value
  @param akwdset - pointer to class which contains set of aliases keywords and corresponding integer values

  @retval 0 - enum integer was not found
  @retval 1 - enum integer was successfully scanned

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long matchekwdint(int *apar, kwdset *akwdset)
{
  long i;
  long ret=0;
  for(i=0; i<akwdset->n; i++)
  {
    if (*apar == akwdset->set[i].id)
    {
      ret = 1;
      break;
    }
  }
  return ret;
}



/**
  Function searches for the given composed integer value apar in alias integers. (used for enum conversion).
  It is supposed that particular aliases in the enum definition akwdset are just power of 2.
  
  Parameters:
  @param apar - searched integer value
  @param akwdset - pointer to class which contains set of aliases keywords and corresponding integer values

  @retval 0 - enum integer was not found
  @retval 1 - enum integer was successfully scanned

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long matchoptekwdint(int *apar, kwdset *akwdset)
{
  long i;
  long ret=1L;
  int aux=*apar;

  for(i=0L; i<akwdset->n; i++)
    aux &= ~(akwdset->set[i].id); // switch off bits for particular aliases

  if (aux) // some remainings in the bit representation
    ret = 0L;

  return ret;
}



/**
  Function checks for errors during enum conversion.
  

  Parameters:
  @param in - pointer to the structure XFILE with opened text file for reading
  @param ret - returned value from the scanf function
  @param reqret - required returned value from the scanf function which was excepted
  @param war - warning indicator (0=warning ignored/ 1=warning / 2=warning considered as error)
  @param akwdset - pointer to class which contains set of aliases keywords and corresponding integer values

  Following parameters are used only for error message 
  @param fmt - format string
  @param idt - index of actually processed conversion in fmt
  @param ptype   - pointer to the actual conversion
  @param subfmt  - string for subformat
  @param errfile - string with source file name where the error was detected
  @param errline - line number in source file where the error was detected
  @param errfunc - string with function name where the error was detected

  @retval 0 - enum integer was not found
  @retval 1 - enum integer was successfully scanned

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long checkenumerr(XFILE *in, long ret, long reqret, long war, kwdset *akwdset, 
                  const char *fmt, char *subfmt, long idt, const char *ptype, const char *errfile,
                  int errline, const char *errfunc)
{
  long i;
  const char *partmsg;

  if (ret != reqret)
  {
    switch (ret)
    {
      case 2:
        partmsg = "cannot read integer id of enum";
        break;
      case 3:
        partmsg = "enum integer id was not found in keyword set";
        break;
      case 4:
        partmsg = "cannot read enum keyword nor integer id";
        break;
      case 5:
        partmsg = "keyword set was not specified in source code file";
        break;
      default:
        partmsg = "cannot read enum";
    }
    getsubfmt(ptype, subfmt, 0);
    print_err("%s\n"
              "(format '%s', subformat '%s', parameter id=%ld)\n"
              "input file: %s, line=%ld, column=%ld", errfile, errline, errfunc, 
              partmsg, fmt, subfmt, idt+1, in->fname, in->line, in->col);
    if (akwdset == NULL)
      fprintf(stderr, "No keyword set was specified\n");
    else
    {
      fprintf(stderr, "Following values are allowed for this enum:\n");       
      for (i=0; i < akwdset->n; i++)
        fprintf(stderr, "%s (%d)\n", akwdset->set[i].alias, akwdset->set[i].id);
    }
    return 1;
  }
  if (war && in->warning)
  {
    getsubfmt(ptype, subfmt, 0);
    if (in->warning == 2)
    {   
      print_err("\n(format '%s', subformat '%s', parameter id=%ld)\n"
                "input file: %s, line=%ld, column=%ld", errfile, errline, errfunc, 
                fmt, subfmt, idt+1, in->fname, in->line, in->col);
      return 1;
    }
    if (in->warning == 1)
    {   
      print_warning("unable to read integer\n(format '%s', subformat '%s', parameter id=%ld)\n"
                    "input file: %s, line=%ld, column=%ld", errfile, errline, errfunc, 
                    fmt, subfmt, idt+1, in->fname, in->line, in->col);
    }
  }
  return 0;
}

/**
  Function performs case sensistive/insensitive searching of substring 'needle' in the string 'haystack', 
  The terminating '\0' characters are not compared. The case sensitivity is controlled by the parameter
  ignorecase: 0 = case sensitive / 1 = case insensitive.

  Parameters:
  @param haystack - searched in string
  @param needle - searched for string
  @param ignorecase - flag for case insensitivity (=1)/ case sensitivity (=0)

  @return The function returns a pointer to the beginning of the substring, or NULL if the substring is not found

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
char *strstr(char *haystack, const char *needle, long ignorecase)
{
  char *ret;
  if (ignorecase==1)
    ret = strstrcis(haystack, needle);
  else
    ret = strstr(haystack, needle);
  return ret;
}



/**
  Function performs case sensistive/insensitive searching of substring 'needle' in the string 'haystack', 
  The terminating '\0' characters are not compared. The case sensitivity is controlled by the parameter
  ignorecase: 0 = case sensitive / 1 = case insensitive.

  Parameters:
  @param haystack - searched in string literal
  @param needle - searched for string
  @param ignorecase - flag for case insensitivity (=1)/ case sensitivity (=0)

  @return The function returns a pointer to the beginning of the substring, or NULL if the substring is not found

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
const char *strstr(const char *haystack, const char *needle, long ignorecase)
{
  const char *ret;
  if (ignorecase==1)
    ret = strstrcis(haystack, needle);
  else
    ret = strstr(haystack, needle);
  return ret;
}



/**
  Function performs case insensitive searching of substring 'needle' in the string 'haystack', The terminating
  '\0' characters are not compared.

  Parameters:
  @param haystack - searched in string
  @param needle - searched for string

  @return The function returns a pointer to the begning of the substring, or NULL if the substring is not found

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
char *strstrcis(char *haystack, const char *needle)
{
  size_t i, j, lh, ln;
  lh = strlen(haystack);
  ln = strlen(needle);
  j = 0;

  if ((ln == 0) || (needle == NULL) || (haystack == NULL))
    return NULL;

  for (i = 0; i<lh; i++)
  {
    if (tolower(haystack[i]) == tolower(needle[j]))
    {
      j++;
      if (j == ln)
        return (char *)(haystack+i-(j-1));
    }
    else
      j = 0;  
  }

  return NULL;
}



/**
  Function performs case insensitive searching of substring 'needle' in the string 'haystack', The terminating
  '\0' characters are not compared.

  Parameters:
  @param haystack - searched in string literal
  @param needle - searched for string

  @return The function returns a pointer to the begning of the substring, or NULL if the substring is not found

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
const char *strstrcis(const char *haystack, const char *needle)
{
  size_t i, j, lh, ln;
  lh = strlen(haystack);
  ln = strlen(needle);
  j = 0;

  if ((ln == 0) || (needle == NULL) || (haystack == NULL))
    return NULL;

  for (i = 0; i<lh; i++)
  {
    if (tolower(haystack[i]) == tolower(needle[j]))
    {
      j++;
      if (j == ln)
        return (haystack+i-(j-1));
    }
    else
      j = 0;  
  }

  return NULL;
}



/**
  Function skips rest of actual line in file in.

  Parameters:
  @param in - pointer to the structure XFILE with opened text file for reading

  @return The function does not return anything.

  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
  modified  11.9.2013,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void skipline(XFILE *in)
{
  in->col = long(strlen(in->lnstr)+1);
  return;
}



/**
  Function reads the ordinary format characters contained in the format substring fmt from the file in.
  Ordinary format characters are defined as characters which does not belong to the
  conversion specification. 
  For example, the format string  "z= %le, ndof= %ld " has to be
  divided into three format substrings  
   -  s1 = "z= ", 
   -  s2 = "%le, ndof= " and 
   -  s3 = "%ld ".

  Any of s1, s2 or s3 can be passed as the argument fmt and fmt_end has to be pointer to the next substring
  i.e. s2, s3 or NULL.

  Following ordinary format substrings will be detected: 
   -  "z= "       for s1,
   -  ", ndof= "  for s2 and
   -  " "         for s3.

  @param in      - pointer to the structure XFILE with opened text file for reading
  @param fmt     - format substring containing at most one conversion specifier 
  @param fmt_end - pointer to the end of format substring fmt or NULL if the fmt is terminated by '\0'
  @param lnstr   - line buffer of the file in
  @param aptr    - pointer of the actual position in lnstr

  @retval NULL - if the ordinary format has been read successfully
  @retval afmt - if a matching failure has been detected (afmt is string with format which caused the failure)

  Created 04.2011,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz  
*/
const char *proc_ord_fmt(XFILE *in, const char *fmt, const char *fmt_end, char *lnstr, char *&astr)
{
  long l, err, ret;
  unsigned int br;
  char *subfmt; // pointer to the one conversion format including surrounding ordinary format characters
                // the first character of the format has to be '%' otherwise no conversion has to be contained in the subfmt
  const char *afmt;   // pointer to the actual position of the ordinary format in the subfmt
  const char *caux;
  char *aux;  


  if (fmt_end)
    l = long(fmt_end - fmt);
  else
    l = long(strlen(fmt)); 
  subfmt = new char[l+1];
  strncpy(subfmt, fmt, l);
  subfmt[l] = '\0';

  if (subfmt[0] == '%')
  {
    afmt = caux = strpbrk(subfmt, "acCdeEfFgGikmnopsSuxX[");  // check for known conversions
    if (*afmt == '[') // for "%[...]" conversion the terminating character ']' has to be found
    {                
      caux = strchr(afmt+1, ']');
      if (afmt[1] == ']')  // if the ']' should be included in the set it has to be the first character after opening '['
        caux = strchr(afmt+2, ']');  // scanning for the terminating ']' has to be started from the third character
      if (afmt[1] == '^')  // if the ']' should be excluded from the set it has to be the first character after '^'
        caux = strchr(afmt+3, ']');  // scanning for the terminating ']' has to be started from the fourth character
    }
    if (caux == NULL)
    {
      print_err("invalid format %s - cannot find terminating ]", __FILE__, __LINE__, __func__, subfmt);
      delete [] subfmt;
      return NULL;
    }
      
    afmt = caux+1; // pointer to the first character after the conversion type specifier (i.e. the first character of the ordinary format)
  }
  else
   afmt = subfmt; // pointer to the ordinary format which begins at the start character of the format string


  l = long(strlen(afmt));
  if (l == 0)  // no ordinary format characters continues after format type
  {
    delete [] subfmt;
    return NULL;
  }
  aux = new char[l+3]; // additional space for %n and '\0' is needed
  do
  {
    if (isspace(*afmt))         // one or more blanks are in the actual format position
    {
      while (isspace(*afmt))    // skip blanks in the format
        afmt++; 

      // skip blanks in the file
      ret = skip_space(in, lnstr, astr);

      if (ret == 2)               // end of file/section was detected and warnings are considered as an errors
      {
        afmt = fmt + (afmt - subfmt);
        delete [] subfmt;
        delete [] aux;
        return afmt;
      }

      if (ret == 1)
      { 
        if (*afmt)  // end of file/section and something left to be read in the ordinary format
        {
          afmt = fmt + (afmt - subfmt);
          delete [] subfmt;
          delete [] aux;
          return afmt;                   
        }
        else                        // end of file/sectio and nothing else to be read
        {
          delete [] subfmt;
          delete [] aux;
          return NULL;
        }
      }
    }
    if (sscanf(afmt, "%s", aux)==1)  // sequence of nonblank oridinary format characters detected in the actual format string
    { 
      l = long(strlen(aux));
      sprintf(aux+l, "%%n");
      br = 0;
      // try to read non-blanks
      err = sscanf(astr, aux, &br);
      if ((err != EOF) && (br == l))  // non-blanks have been read successfully 
      {
        afmt += l;     // change the actual position in fmt
        astr += l;     // change the position in lnstr
        in->col += l;  // change the column number
      }
      else  // matching failure has been detected 
      {        
        afmt = fmt + (afmt - subfmt);
        delete [] subfmt;
        delete [] aux;
        return afmt;
      }
    }
  } while(*afmt != '\0');

  delete [] subfmt;
  delete [] aux;
  return NULL;
}



/**
  Function skips the blank characters ' ', '\n', '\r', \t', '\v' 
  (defined in the standard isspace function) in the file in.

  @param in    - pointer to the structure XFILE with opened text file for reading
  @param lnstr - line buffer of the file in
  @param astr  - pointer of the actual position in lnstr
  
  @retval 0 - on success
  @retval 1 - on end of file/section
  @retval 2 - on maximum line length error

  Created 04.2011,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz  
*/
long skip_space(XFILE *in, char *lnstr, char *&astr)
{
  long l;
  int  br;
  long pos;
  int  fbr;
  char rfmt[20];   // format string for reading of new line using fscanf

  do
  {
    l  = long(strlen(astr));
    br = 0;
    if (l != 0)
    {
      sscanf(astr, " %n", &br); // scan ' ', '\t', '\v' from the actual position of line
      astr += br;               // change the actual postion in the line by the number of bytes read in the above sscanf 
      in->col += br;            // change the column number by the number of bytes read in the above sscanf 
      if (br != l)              // if the length of astr does not match to bytes read ->
        return 0;               // nonblank characters were matched in the astr
    }

    // read new line character '\n'
    // ----------------------------
    // check end of file
    if (feof(in->file))
      return 1;
    // check end of section
    pos = ftell(in->file);
    if (in->index_created)
    {
      if ((in->line == in->asect->end_secln) && 
          (pos > in->asect->end_secpos))
        return 1;
    }
    fscanf(in->file, "%*1[\n]"); // read terminating \n
    in->line++;
    in->col=1;
    in->lnfpos = ftell(in->file);// backup of beginning line position

    // read new line
    // -------------
    memset(lnstr, 0, sizeof(*lnstr)*in->give_maxlnsize()+1);
    // check end of file
    if (feof(in->file))
      return 1;

    // check end of section
    pos = ftell(in->file);
    if (in->index_created)
    {
      if ((in->line == in->asect->end_secln) && 
          (pos > in->asect->end_secpos))
        return 1;
    }

    // prepare read format string for fscanf function in the form "%x[^\n]", where x is in->maxlnsize
    sprintf(rfmt, "%%%ld[^\n]%%n", in->give_maxlnsize());
    fbr = 0;
    fscanf(in->file, rfmt, lnstr, &fbr); //read new line without terminating \n
    if (check_maxlnsize(in, fbr))
      return 2;

    // cut string before end of actual section
    cut_str_sec(in, lnstr);
    // removing trailing \r form the last character of string (due to DOS \r\n newline)
    l = long(strlen(lnstr));
    if (l > 0)
    {
      if (lnstr[l-1] == '\r')
        lnstr[l-1] = '\0';
    }
    clearcomments(lnstr);   // clear possible comments
    astr = lnstr;       // setup of actual position in lnstr
  } while(1);

  return 0;    
}



/** 
  The function checks whether the '\n' character is included in the 
  accepted set of characters for format "%[...]".

  @param fmt - sreached format string

  @retval 0 - the '\n' character is not included in the format set
  @retval 1 - the '\n' character is included in the format set

  Created 04.2011 by Tomas Koudelka
*/
long check_newline_fmt(const char *fmt)
{
  long caretflag = 0;
  const char *aux = strchr(fmt, '[');
  aux += 1;

  if (*aux == '^')
  {
    caretflag = 1;  // excluded characters are specified in the format "%[...]"
    aux++;
  }
  
  while (*aux)
  {
    if ((*aux == '\n') && (caretflag == 0))// backslash character is encountered
      return 1;
    if ((*aux == '\n') && (caretflag == 1))// backslash character is encountered
      return 0;
    aux++;  // process next character in the string aux
  }
  if (caretflag) // if no \n was found in the excluded character list -> \n is accepted by the format "%[...]"
    return 1;
  return 0;
}



/**
  The function searches an occurence of the connecting character '@' in the string astr. 
  The connecting character is searched from the end of astr and only blank characters may
  precede it. It is used for the detection of multiline strings. If the connecting '@' is found then
  it is replaced by '\n'
  
  @param astr - pointer to the string terminated by \0
  
  @retval 0 - no connecting character was found at the end of string
  @retval 1 - connecting character was found at the end of string

  Created by Tomas Koudelka, 08.2011   
*/
long detect_multln(char *astr)
{
  long i, l;
  l = long(strlen(astr));
  for(i=l-1; i >= 0; i--)
  {
    if (isspace(astr[i]))   continue;
    if (astr[i] == '@')
    {
      astr[i] = '\n'; 
      return 1;
    }
  }  
  return 0;
}




/**
  The function checks whether the %a format skips the initial white space characters including \r and \n.
  If the format looks like %a , %#a , %*a or %*#a where # is an integer number than it matches
  even empty strings. If the format type contains whitespaces between '%' and conversion 
  specifier (% a , %* a, % #a or %* #a) only non-empty strings can be matched and the initial 
  whitespaces will be skipped.
  

  @param fmt - format string starting with %a conversion

  @retval 0 - the %a format matches even empty strings
  @retval 1 - the %a format does not match empty strings

  Created by Tomas Koudelka, 03.2012
*/
long a_fmt_skips_whitespaces(const char *fmt)
{
  if (fmt[0] != '%')
  {
    print_err("format string does not start with '%' character", __FILE__, __LINE__, __func__);
    abort();
  }
  if (fmt[1] == '*')
  {
    if (fmt[2] == ' ')
      return 1;
    else
      return 0;
  }
  if (fmt[1] == ' ')
    return 1;
  else
    return 0;

  return 0;
}



/**
  Function decomposes name of file(with path) to path,filename and suffix.
  ppp/ppp/fff.fff.sss -> path = "ppp/ppp" , name = "fff.fff" , suffix = ".sss"
   
  @param file - decomposed name of file
  @param path,name,suffix - components of name of file (output)
  
  @return The function returns pointers to decomposed parts in the parameters
          path, name and suffix for which the memory was allocated.
 
  Created  by Ladislav Svoboda, 1.4.2003, termit@cml.fsv.cvut.cz
 */
void filename_decomposition (const char *file,char *&path,char *&name,char *&suffix)
{
  char *fix;
  path = new char[strlen(file)+1];
  name = new char[strlen(file)+1];
  suffix = new char[strlen(file)+1];
  
  strcpy (path,file);
  fix = strrchr(path,'/');
  if (fix==NULL){
    strcpy (name,path);
    path[0] = '\0';
  }
  else{
    strcpy (name,fix+1);
    fix[1] = '\0';
  }
  
  fix = strrchr(name,'.');
  if (fix==NULL)
    suffix[0] = '\0';
  else{
    strcpy (suffix,fix);
    fix[0] = '\0';
  }
}



/**
  Function decomposes name of file(with path) to path, filename and extension.
  r:/ppp/ppp/fff.fff.sss -> path = "r:/ppp/ppp" , name = "fff.fff" , ext = ".sss"
   
  @param[in] fnstr - pointer to the string of decomposed file name
  @param[out] path - pointer into the fnstr to the path beginning or NULL if no path is given
  @param[out] name - pointer into the fnstr to the file name beginning
  @param[out] suffix - pointer into the fnstr to the file extension beginning or NULL if no extension is given
  
  @return The function returns pointers to decomposed parts in the parameters
          path, name and ext, no memory is allocated for these pointers.
 
  Created by Tomas Koudelka, 10.2023
 */
void filename_decomp_ptr (const char *fnstr, const char *&path, const char *&name, const char *&ext)
{
  path = strrchr(fnstr,'/');
  if (path == NULL)
    path = strrchr(fnstr,'\\');

  if (path)
    name = path + 1;
  else
    name = fnstr;

  ext = strrchr(name, '.');
}



/**
  Function decomposes name of file(with path) to path, filename and extension.
  r:/ppp/ppp/fff.fff.sss -> path = "r:/ppp/ppp" , name = "fff.fff" , ext = ".sss"
   
  @param[in] fnstr - pointer to the string of decomposed file name
  @param[out] path - pointer into the fnstr to the path beginning or NULL if no path is given
  @param[out] name - pointer into the fnstr to the file name beginning
  @param[out] suffix - pointer into the fnstr to the file extension beginning or NULL if no extension is given
  
  @return The function returns pointers to decomposed parts in the parameters
          path, name and ext, no memory is allocated for these pointers.
 
  Created by Tomas Koudelka, 10.2023
 */
void filename_decomp_ptr (char *fnstr, char *&path, char *&name, char *&ext)
{
  path = strrchr(fnstr,'/');
  if (path == NULL)
    path = strrchr(fnstr,'\\');

  if (path)
    name = path + 1;
  else
    name = fnstr;

  ext = strrchr(name, '.');
}



/**
  The function extracts path from the filename given by parameter file.

  @param file[in] - string pointer with filename
  @param p[in/out] - pointer to preallocated string where the returned path will be stored in. If the file argument contains 
                     just the file name without any path, the resulting string will contain actual directory path "./".
                     If the pointer is NULL on input, then new memory is being allocated in the function and the caller 
                     is responsible for its deallocation with the help of delete [] operator.

  @return The function returns pointer to the string with path, i.e. value of argument p.

  Created by Tomas Koudelka, 08/2021.
*/
char* get_path (char *file, char *&p)
{
  size_t l = 0;
  char *fix;
  if (p == NULL)
    p = new char[strlen(file)+3];
  
  fix = strrchr(p, '/');
  if (fix == NULL){
    fix = strrchr(file, '\\');
  }
  if (fix == NULL){
    l = 2;
    strncpy(p, "./", l);
  }    
  else{
    l = fix+1-file;
    strncpy(p, fix, l);
  }
  fix[l] = '\0';
  return p;
}
