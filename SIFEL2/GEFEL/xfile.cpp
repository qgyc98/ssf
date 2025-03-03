#include <string.h>
#include "xfile.h"



XFILE::XFILE()
{ 
  file = NULL; 
  line = col = lnfpos = 0L;
  fname = lnstr = NULL;

  warning = 0; 
  maxlnover = 0L;

  kwdmode = ignore_kwd;
  ignorecase = 0L;
  
  index_created = id_sec = 0L; 
  num_sec = 0L;
  asect = sect = NULL;
 
  lnpref = lnpostf = NULL;
  lnprefl = lnpostfl = 0L;
  
  maxlnsize = 0;
}



XFILE::~XFILE()
{
  delete [] lnstr;
  delete [] sect;
  delete [] fname;
  if (file) fclose(file);
}


/**
  The function returns actual value of maximum line length i.e. size 
  of line buffer for the given XFILE.

  @return Returns the maximum accepted line length.

  Created by Tomas Koudelka, 11.9.2013
*/
long XFILE::give_maxlnsize() 
{
  return maxlnsize;
}



/**
  The function sets new value of maximum line length i.e. size 
  of line buffer for the given XFILE. The size of line buffer is given by 
  parameter lnsize.

  @param lnsize - size of line buffer

  @return The function does not return anything but it changes 
          maxlnsize and lnstr data members.

  Created by Tomas Koudelka, 11.9.2013
*/
void XFILE::set_maxlnsize(long lnsize)
{
  maxlnsize = lnsize;
  if (lnstr)
    delete [] lnstr;
  lnstr = new char[maxlnsize+1];
  memset(lnstr, 0, sizeof(*lnstr)*maxlnsize);
  return;
}



/**
  The function sets the position in file to the end of actual line but 
  before the the \n. In the next call of xfscanf, a new line will have to be read. 

  @return

  Created by Tomas Koudelka, 29.4.2016
*/
void XFILE::set_eoln()
{
  col = long(strlen(lnstr)+1);
}
