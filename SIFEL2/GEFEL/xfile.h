#ifndef XFILE_H
#define XFILE_H
#include <stdio.h>
#include "kwdset.h"

/** aliases for keyword handling - used in XFILE structure and xfscanf function:
    ignore_kwd       - %k conversions are ignored, no keywords are searched for/printed out
    sequnet_mode     - keywords are searched sequentially as the XFILE is being processed, 
                       the searching starts from the actual position in the XFILE and the next string 
                       must match to the given keyword otherwise the error code is returned
    line_mode        - keywords are searched within the scope of the acutal line, if the keyword is not found
                       the error code is returned
    sect_mode_seq    - keywords are searched within the scope of the actual section starting from the actual 
                       position in the section and the next string must match to the given keyword 
                       otherwise the error code is returned. It is similar mode as sequent_mode but the ends
                       of sections are checked for.
    sect_mode_full   - keywords are searched within the scope of the actual section starting from the beginning of 
                       the actual section for each %k conversion in the format string. If the given keyword was not found
                       in the given section the error code is returned.
    sect_mode_ignore - same as the ignore mode but ends of section are checked for
    sect_mode_fwd    - keywords are searched within the scope of the actual section starting from the actual 
                       position in the section. The keywords are searched for until the end of section is reached.
                       If the keyword was not found from the actual position in the section till the end of section 
                       then the error code is returned
*/
enum kwd_handling {ignore_kwd=0, sequent_mode = 1, line_mode=2, sect_mode_seq = 3, sect_mode_full = 4, sect_mode_ignore = 5, sect_mode_fwd = 6};



struct xfsection;
struct XFILE;



/**
  Structure of file description. This file type is used for xfscanf function and 
  it has extended functionality comparing with traditional FILE type - 
   - actual line and column number, 
   - file name string, 
   - error and warning management,
   - support for keywords.

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz,  12.2006
*/
struct XFILE
{
  FILE *file;     ///< pointer to standard file structure
  long line;      ///< actual line
  long col;       ///< actual column
  long lnfpos;    ///< actual position of actual line beginning in file
  char *fname;    ///< name of opened file
  char *lnstr;    ///< line buffer where the actual file line is stored


  int warning;    ///< flag which handles warning (0=no warning/1=print warning/2=treat warning as an error)
  long maxlnover; ///< flag with the number of last oversized line (=0 if no oversized line detected)

  kwd_handling kwdmode; ///< keyword handling mode  
  long ignorecase;      ///< flag for case sensitive(=0)/case insensitive(=1) keyword searching 

  // following data members are necessary only in case that kwdmode=sect_mode
  long index_created; ///< flag for created index of particular sections
  long id_sec;        ///< index of actual section
  xfsection *asect;   ///< pointer to actual section
  long num_sec;       ///< number of detected sections
  xfsection *sect;    ///< array with descriptions of detected sections  

  char *lnpref;       ///< prefix string of line beginning (for XFILE output)
  char *lnpostf;      ///< postfix string of line end (for XFILE output)
  long lnprefl;       ///< length of lnpref string
  long lnpostfl;      ///< length of lnpostf string

  private:
   long maxlnsize;    ///< max line size
  public:
   XFILE();
  ~XFILE();
   long give_maxlnsize();           ///< returns size of line buffer
   void set_maxlnsize(long lnsize); ///< sets size of line buffer
   void set_eoln();                 ///< sets actual position in line buffer to the end
};



/**
  Structure contains data about XFILE sections. It is used 
  in case that kwdmode=sect_mode

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz,  12.2006
*/
struct xfsection
{
  long beg_secpos;    ///< position of section beginning in xfile (section begins at the first character after the keyword for the section begining)
  long end_secpos;    ///< position of section end in xfile  (section ends at the last character befor the keyword for the section end)
  long beg_secln;     ///< line number of section beginning
  long end_secln;     ///< line number of section end
  long beg_seccol;    ///< column number of section begining (immediately after keyword of section beginning)
  long end_seccol;    ///< column number of section end (the last position before keyword of section end)
  enumstr name;       ///< section name and enum alias
};
#endif
