#include "intools.h"
#include "iotools.h"
#include "xfile.h"
#include <string.h>
#include <ctype.h>
/**
 This file contains functions for processing file with keywords and comments blocks
*/

static long Ngkw = 2; ///< number of general keywords
static const char *Genkw[2] = {"begin", "end"}; ///< general keywords string array
static long Npkw = 21; ///< number of property keywords
static const char *Propkw[21] = {"ndofn",    "bocon", "nload",  "crsec",   "lcsys",  "eltype",
                           "mater",    "eload", "dloadn", "dloadel", "inicon", "comcn",
                           "loadedge", "sscomp", "ntemp", "tfunct", "loadsurf", "elsurfload", 
                           "ntimefunc", "etimefunc", "bocsurf"};
                           ///< property keywords string array

/**
  Function comparess string s with Genkw array and returns result

  @param s   - pointer to compared string.
  @param ptr - output parameter which contains pointer to array Genkw where
               is same keyword as in the s string or NULL when no general
               keyword matches.

  @retval -1 - no general keyword matches string s.
  @retval  i - i-th general keyword matches string s. Indexing starts from 0

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long isgkw(char *s, char *&ptr)
{
  ptr = NULL;
  for (long i = 0; i < Ngkw; i++)
  {
    if ((ptr = strstr(s, Genkw[i])) != NULL)
      return(i);
  }
  return(-1);
}

/**
  Function comparess string s with Propkw array and returns result

  @param s   - pointer to compared string.
  @param ptr - output parameter which contains pointer to array Propkw where
               is same keyword as in the s string or NULL when no property
               keyword matches.

  @retval -1 - no property keyword matches string s.
  @retval  i - i-th property keyword matches string s. Indexing starts from 0

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long ispkw(char *s, char *&ptr)
{
  for (long i = 0; i < Npkw; i++)
  {
    if ((ptr = strstr(s, Propkw[i])) != NULL)
      return(i);
  }
  return(-1);
}

/**
  Function searches string s for comment character # and returns result

  @param s   - pointer to searched string.
  @param ptr - output parameter which contains pointer to string s where
               was found character # or NULL when no # chracter found

  @retval 1 - There was found # character in the string.
  @retval 0 - No # chracterwas was found.

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long iscomment(char *s, char *&ptr)
{
  if ((ptr = strstr(s, "#")) != NULL)
    return(1);
  return(0);
}

/**
  Function deletes part of string s behind comment character # and puts \0 chracter
  instead it. If only blank chracters precede #, then \0 is placed to the first (0-th)
  position of string

  @param s   - pointer to searched string.

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
void delcomment(char *s)
{
  long l = long(strlen(s));
  long  nonsp = -1;

  for (long i = 0; i < l; i++)
  {
    if (! isspace(s[i]))
      nonsp++;
    if (s[i] == '#')
    {
      if (nonsp)
        s[i] = '\0';
      else
        s[0] = '\0';
      break;
    }
  }
}

/**
  Function deletes comment blocks from file one block denoted by general keywords 'begin' and 'end'.
  This block is rewritten without comments to the temporary file.

  @param in - pointer to opened file with block which will be cleaned out of comments

  @retval pointer to opened temporary file where is the result of cleaning process.

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
FILE *cleancommblock(FILE *in)
{
  FILE *ret;
  long posb, curpos, comfl;
  size_t pose;
  int c;

  ret = tmpfile();
//  ret = fopen("temppd","w+");
  if (ret == NULL)
    return(ret);
  getgkwid(in, 0);  // searches begin keyword
  getnextln(in);    // moves on the next line, where begin commented problem description
  posb = ftell(in); // remeber position of block's end
  getgkwid(in, 1);  // searches end keyword
  pose = ftell(in); // remeber position of block's end
  pose -= strlen(Genkw[1]); //position of the first character 'end' keyword
  fseek(in, posb, SEEK_SET); // psoition of the first charcter of problem description block
  curpos = posb;
  comfl = 0;
  while (size_t(curpos) < pose)
  {
    c = fgetc(in);
    if (c == '#')
    {
      getnextln(in);
      comfl = 1; // after skipping comment lines a '\n' must be written to the ret
    }
    else
    {
      if (comfl)
        fputc('\n', ret);
      fputc(c, ret);
      comfl = 0; // not neccesarry write '\n' due to comment lines
    }
    curpos = ftell(in);
  }
  rewind(ret); // begin of file
  return(ret);
}

/**
  Function searches file given by parameter in for \n character and sets file position
  indicator on it.

  @param in - pointer to opened file where \n is searched

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
void geteoln(FILE *in)
{
  int c = 0;
  while ((! feof(in)) && ((c = fgetc(in)) != '\n'));
  if (c == '\n')
    fseek(in, -1, SEEK_CUR);
}

/**
  Function searches file given by parameter in for \n character and sets file position
  indicator after it.

  @param in - pointer to opened file where \n is searched

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
void getnextln(FILE *in)
{
  int c = 0;
  while ((! feof(in)) && ((c = fgetc(in)) != '\n'));
}

/**
  Function searches file given by parameter in from current position for WHITESPACE characters
  and sets file position indicator to the first nonwhitespace character i.e. it searches next text.

  @param in - pointer to opened file where \n is searched

  created  9.2003,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 9.2003, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
void getnexttxt(FILE *in)
{
  int c = 0;
  while (! feof(in))
    {
      c = fgetc(in);
      if (! isspace(c))
	{
	  //fseek(in, SEEK_CUR, -1);
	  ungetc(c, in);
	  break;
	}
    };
}

/**
  Function searches file given by parameter in for id-th general keyword.
  Comment blocks are accepted and skipped. On success, the file position indicator is
  set after last character searched keyword.

  @param in - pointer to opened file where general keyword is searched
  @param id - index of desired general keyword

  @retval 1 - if fails reading of input file
  @retval 0 - on succes

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long getgkwid(FILE *in, long id)
{
  char buf[1025]; char *ptr, *ptrc;
  long pos = ftell(in);
  long idkw;
  int  cb;
  do
  {
    memset(buf, 0, sizeof(*buf)*1025);
    if (fscanf(in, "%1024s%n", buf, &cb) != 1)
    {
      fprintf(stderr, "\nError searching general keyword ""%s"" in function getgkwid()", Genkw[id]);
      fprintf(stderr, "\n in file input.cpp\n");
      return(1);
    }
    pos += cb;
    idkw = isgkw(buf, ptr);
    // search eventually comment
    if (iscomment(buf, ptrc))
    {
      if ((ptrc < ptr) || (idkw < 0))
      {
        pos -= cb; // position before last buf reading
        pos += long(ptrc - buf); // position before the first buf character
        pos += long(strlen("#")); // position of last character of the desired keyword
        pos++; // position after last character of desired keyword
        fseek(in, pos, SEEK_SET);
        geteoln(in);
        pos = ftell(in);
        idkw = -1;
      }
    }
  }
  while (idkw != id);
  pos -= cb; // position before last buf reading
  pos += long(ptr - buf); // position before the first buf character
  pos += long(strlen(Genkw[id])); // position of last character of the desired keyword
  pos++; // position after last character of desired keyword
  fseek(in, pos, SEEK_SET);
  return(0);
}


/**
  Function searches file given by parameter in for id-th general keyword.
  Comment blocks are accepted and skipped. On success, the file position indicator is
  set after last character searched keyword.

  @param in - pointer to opened file where general keyword is searched
  @param id - index of desired general keyword
  @param numl - number of processed lines

  @retval 1 - if fails reading of input file
  @retval 0 - on succes

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long getgkwid(XFILE *in, long id, long &numl)
{
  char *buf; char *ptr;
  long pos;
  long idkw;
  unsigned int  cb;

  buf = new char[in->give_maxlnsize()];
  memset(buf, 0, sizeof(*buf)*in->give_maxlnsize());
  numl = -1;
  do
  {
    in->lnfpos = pos = ftell(in->file);
    if ((cb = inputln(in->file, buf, in->give_maxlnsize())) == 0)
    {
      fprintf(stderr, "\nError searching property keyword in function getpkw()");
      fprintf(stderr, "\n in file input.cpp\n");
      numl = 0;
      return(1);
    }
    numl++;
    clearcomments(buf);
    idkw = isgkw(buf, ptr);
  }
  while (idkw != id);
  pos += long(ptr - buf); // position before the first buf character
  pos += long(strlen(Propkw[idkw])); // position of last character of the desired keyword
  fseek(in->file, pos, SEEK_SET);
  in->col = pos - in->lnfpos+1;
  delete [] buf;
  return(0);
}

/**
  Function searches file given by parameter in for id-th property keyword.
  Comment blocks are accepted and skipped. On success, the file position indicator is
  set after last character searched keyword.

  @param in - pointer to opened file where property keyword is searched
  @param id - index of desired property keyword

  @retval 1 - if fails reading of input file
  @retval 0 - on succes

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long getpkw(FILE *in)
{
  char buf[1025], *ptr, *ptrc;
  long pos = ftell(in);
  int  cb;
  long idkw;
  do
  {
    if (fscanf(in, "%1024s%n", buf, &cb) != 1)
    {
      fprintf(stderr, "\nError searching property keyword in function getpkw()");
      fprintf(stderr, "\n in file input.cpp\n");
      return(1);
    }
    pos += cb;
    idkw = ispkw(buf, ptr);
    // search eventually comment
    if (iscomment(buf, ptrc))
    {
      if ((ptrc < ptr) || (idkw < 0))
      {
        pos -= cb; // position before last buf reading
        pos += long(ptrc - buf); // position before the first buf character
        pos += long(strlen("#")); // position of last character of the desired keyword
        pos++; // position after last character of desired keyword
        fseek(in, pos, SEEK_SET);
        geteoln(in);
        pos = ftell(in);
        idkw = -1;
      }
    }
  }
  while (idkw < 0);
  pos -= cb; // position before last buf reading
  pos += long(ptr - buf); // position before the first buf character
  pos += long(strlen(Propkw[idkw])); // position of last character of the desired keyword
  pos++; // position after last character of desired keyword
  fseek(in, pos, SEEK_SET);
  return(idkw);
}



/**
  Function searches file given by parameter in for property keywords and counts its number.
  File is proccessing until any general keyword found.
  Comment blocks are accepted and skipped. After searching the file position indicator is
  set to the original position before searching.

  @param in - pointer to opened file where property keywords are searched
  @param numl - number of lines found during the searching

  @retval 1 - if fails reading of input file
  @retval 0 - on succes

  created  2.2007,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long getnprop(XFILE *in, long &numl)
{
  long ret = 0;
  char *buf, *ptr, *ptrg;
  long origpos = ftell(in->file);
  long idkw, idgkw;
  long cb;

  buf = new char[in->give_maxlnsize()];
  memset(buf, 0, sizeof(*buf)*in->give_maxlnsize());
  numl = -1;
  do
  {
    if((cb = inputln(in->file, buf, in->give_maxlnsize())) == 0)
    {
      fprintf(stderr, "\nError counting property input lines in function getnprop()");
      fprintf(stderr, "\n in file input.cpp\n");
      return(-1);
    }
    numl++;
    clearcomments(buf);
    if ((idkw = ispkw(buf, ptr)) > -1)
      ret++;
    idgkw = isgkw(buf, ptrg);
  }
  while (idgkw < 0);
  fseek(in->file, origpos, SEEK_SET);
  delete [] buf;
  return(ret);
}

/**
  Function searches file given by parameter in for property keywords and counts its number.
  File is proccessing until any general keyword found.
  Comment blocks are accepted and skipped. After searching the file position indicator is
  set to the original position before searching.

  @param in - pointer to opened file where property keywords are searched
  @param id - index of desired property keyword

  @retval 1 - if fails reading of input file
  @retval 0 - on succes

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long getnprop(FILE *in)
{
  long pos, ret = 0;
  char buf[1025], *ptr, *ptrc, *ptrg;
  long origpos = pos = ftell(in);
  long idkw, idgkw;
  int cb;

  do
  {
    if (fscanf(in, "%1024s%n", buf, &cb) != 1)
    {
      fprintf(stderr, "\nError counting property input lines in function getnprop()");
      fprintf(stderr, "\n in file input.cpp\n");
      return(-1);
    }
    pos += cb;
    if ((idkw = ispkw(buf, ptr)) > -1)
      ret++;
    idgkw = isgkw(buf, ptrg);
    // search eventually comment
    if (iscomment(buf, ptrc))
    {
      if ((ptrc < ptr) || (idkw < 0))
      {
        pos -= cb; // position before last buf reading
        pos += long(ptrc - buf); // position before the first buf character
        pos += long(strlen("#")); // position of last character of the desired keyword
        pos++; // position after last character of desired keyword
        fseek(in, pos, SEEK_SET);
        geteoln(in);
        pos = ftell(in);
        if (idkw > -1)
          ret--;
      }
      if ((idgkw > -1) && (ptrc < ptrg))
        idgkw = -1;
    }
  }
  while (idgkw < 0);
  fseek(in, origpos, SEEK_SET);
  return(ret);
}



/**
  Function reads integer number type long from file given by parameter in.
  Comment blocks are accepted and skipped. On success, the file position indicator is
  set after last character searched number.

  @param in  - pointer to opened file where the number is searched
  @param num - output parameter for result

  @retval 1 - if fails reading of input file
  @retval 0 - on succes

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long getlong(FILE *in, long &num)
{
  char buf[1025], *ptrc;
  long pos = ftell(in);
  int  cb, tcb;
  long idl;
  do
  {
    if (fscanf(in, "%1024s%n", buf, &cb) != 1)
    {
      fprintf(stderr, "\nError reading long number in function getlong()");
      fprintf(stderr, "\n in file input.cpp\n");
      return(1);
    }
    pos += cb;
    idl = sscanf(buf, "%ld%n", &num, &tcb);
    // search eventually comment
    if (iscomment(buf, ptrc))
    {
      if ((ptrc < (buf+tcb)) || (idl == 0))
	{
	  /*        pos -= cb; // position before last buf reading*/
	  pos -= long(strlen(buf)); // position of the first buf character
	  pos += long(ptrc - buf); // position of the commentary character
	  fseek(in, pos, SEEK_SET);
	  geteoln(in);
	  pos = ftell(in);
	  idl = 0;
	}
    }
  }
  while (idl == 0);
  /*  pos -= cb; // position before last buf reading
      pos += tcb; // position before the first buf character
      pos++; // position after last character of desired number*/
      pos -= long(strlen(buf)); // position of the first buf character
  pos += tcb; // position after last character of desired number
  fseek(in, pos, SEEK_SET);
  return(0);
}

/**
  Function reads integer number type int from file given by parameter in.
  Comment blocks are accepted and skipped. On success, the file position indicator is
  set after last character searched number.

  @param in  - pointer to opened file where the number is searched
  @param num - output parameter for result

  @retval 1 - if fails reading of input file
  @retval 0 - on succes

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long getint(FILE *in, int &num)
{
  char buf[1025], *ptrc;
  long pos = ftell(in);
  int  cb, tcb;
  long idl;
  do
  {
    if (fscanf(in, "%1024s%n", buf, &cb) != 1)
    {
      fprintf(stderr, "\nError reading int number in function getlong()");
      fprintf(stderr, "\n in file input.cpp\n");
      return(1);
    }
    pos += cb;
    idl = sscanf(buf, "%d%n", &num, &tcb);
    // search eventually comment
    if (iscomment(buf, ptrc))
    {
      if ((ptrc < (buf+tcb)) || (idl == 0))
      {
        pos -= long(strlen(buf));  // position of the first buf character
        pos += long(ptrc - buf); // position of commentary character (#)
        fseek(in, pos, SEEK_SET);
        geteoln(in);
        pos = ftell(in);
        idl = 0;
      }
    }
  }
  while (idl == 0);
  pos -= long(strlen(buf)); // position of the first buf character
  pos += tcb;         // position after last character of desired number
  fseek(in, pos, SEEK_SET);
  return(0);
}

/**
  Function reads real number type double from file given by parameter in.
  Comment blocks are accepted and skipped. On success, the file position indicator is
  set after last character searched number.

  @param in  - pointer to opened file where the number is searched
  @param num - output parameter for result

  @retval 1 - if fails reading of input file
  @retval 0 - on succes

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long getdouble(FILE *in, double &num)
{
  char buf[1025], *ptrc;
  long pos;
  int  cb, tcb;
  long idl;
  // reach the first nonwhitespace character
  do
  {
    pos = ftell(in);
    if (fscanf(in, "%1024s%n", buf, &cb) != 1)
    {
      fprintf(stderr, "\nError reading double number in function getlong()");
      fprintf(stderr, "\n in file input.cpp\n");
      return(1);
    }
    pos += cb;
    idl = sscanf(buf, "%lf%n", &num, &tcb);
    // search eventually comment
    if (iscomment(buf, ptrc))
    {
      if ((ptrc < (buf+tcb)) || (idl == 0))
      {
        pos -= long(strlen(buf));  // position of the first buf character
        pos += long(ptrc - buf); // position of the commentary character
        fseek(in, pos, SEEK_SET);
        geteoln(in);
        pos = ftell(in);
        idl = 0;
      }
    }
  }
  while (idl == 0);
  pos -= long(strlen(buf));  // position of the first buf character
  pos += tcb;          // position after last character of desired number
  fseek(in, pos, SEEK_SET);
  return(0);
}

/**
  Function reads string from file given by parameter in. Initial whitespaces are skipped.
  String str has maximu length maxlen
  chracters including terminating \0.
  Comment blocks are accepted and skipped. On success, the file position indicator is
  set after last character searched string.

  @param in  - pointer to opened file where the string is searched
  @param num - output parameter for result

  @retval 0 always

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long getstring(FILE *in, char *str, long maxlen)
{
  long i = 0;
  int c;

  // reach the first nonwhitespace character
  while (1)
  {
    c = fgetc(in);
    if (feof(in))
      return(1);
    if (! isspace(c))
    {
      ungetc(c, in);
      break;
    }
  }

  // read data to str
  c = 0;
  while (i < maxlen)
  {
    if(feof(in))
      break;
    c = fgetc(in);  
    str[i] = char(c);
    if (c == '\n')   // break if newline chracter has been reached
    {
      ungetc(c, in);
      break;
    }
    if (c == '#') // cutoff comment
    {
      if (i != 0) // if any nonwhitespace chracter was read then break
                  // else continue reading string on the next line
      {
        geteoln(in); // reach end of line
        break;
      }
      else
      {
        getnextln(in);
        continue;
      }
    }
    i++;
  }
  str[i] = '\0'; // terminating chracter
  // reach end of line
  if ((i == maxlen) && (c != '\n'))
    while (! feof(in) && (fgetc(in) != '\n'));
  return(0);
}

/**
  Function reads string from file given by parameter in. String str has maximu length maxlen
  chracters including terminating \0.
  Comment blocks are accepted and skipped. On success, the file position indicator is
  set after last character searched string.

  @param in  - pointer to opened file where the string is searched
  @param num - output parameter for result

  @retval 0 always

  created  12.2001,  Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long getstring2(FILE *in, char *str, long maxlen)
{
  long i = 0;
  int c;

  // read data to str
  c = 0;
  while (i < maxlen)
  {
    if(feof(in))
      break;
    c = fgetc(in);  
    str[i] = char(c);
    if (c == '\n')   // break if newline chracter has been reached
    {
      // ungetc(c, in); newline character have to be read to move to the next line 
      break;
    }
    if (c == '#') // cutoff comment
    {
      if (i != 0) // if any nonwhitespace chracter was read then break
                  // else continue reading string on the next line
      {
        geteoln(in); // reach end of line
        break;
      }
      else
      {
        getnextln(in);
        continue;
      }
    }
    i++;
  }
  str[i] = '\0'; // terminating chracter
  // reach end of line
  if ((i == maxlen) && (c != '\n'))
    while (! feof(in) && (fgetc(in) != '\n'));
  return(0);
}

/**
  Function reads maxlen characters to string anystr from th file in. Function reads
   file till the line treminating character \n and stores chracters to string given by parameter
   anystr until maxlen chractres are read.

 @param in     - pointer to opened file
 @param anystr - output parameter, pointer to allocated string  where read chracters are stored
 @param maxlen - maximu number of read chracters

 @retval number of successfully read bytes
*/
long inputln(FILE* in, char* anystr, long maxlen)
{
  long i = 0;
  char c;

  c = '\0';
  anystr[0] = '\0';
  do
  {
    if (i >= maxlen)
      break;
    if (feof(in))
      break;
    c = char(fgetc(in));
    if ((c == '\n') && strlen(anystr) == 0)
      continue;
    if ((c == '\n') && i != 0)
      break;
    anystr[i] = c;
    i++;
  } while (1);
  anystr[i] = '\0';
  if ((i == maxlen) && (c != '\n'))
    while (! feof(in) && (fgetc(in) != '\n')) ;
  
  return (i);
}
