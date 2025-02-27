#ifndef IOTOOLS_H
#define IOTOOLS_H
#include "xfile.h"
#include "kwdset.h"




/// function prints error message to standard error output device, format specifiers are allowed in emsg
void print_err(const char *emsg, const char *errfile, int errln, const char *errfunc, ...);

/// function prints warning message to standard error output device, format specifiers are allowed in wmsg
void print_warning(const char *wmsg, const char *warnfile, int warnln, const char *warnfunc, ...);

/// function prints error message to standard error output device for parallel executions, format specifiers are allowed in emsg
void par_print_err(int n,const char *emsg, const char *errfile, int errln, const char *errfunc, ...);

/// function prints error message to standard error output device for parallel executions, format specifiers are allowed in emsg
void par_print_err(int n,const char *pname, const char *emsg, const char *errfile, int errln, const char *errfunc, ...);

/// function opens file of XFILE type
XFILE *xfopen(const char *name, const char *mode);

/// function closes file of XFILE type
int xfclose(XFILE *&f);

/// functions searches for keyword in a string for given section in xfscanf function
long getkwd_sect(XFILE *in, char *lnstr, char *&aptr, const char * apar, long detect_allkwd);

/// function loactes keyword in the given string
char *locate_kwd(char *lnstr, const char *kwd, long ignorecase, long &noc);

/// function detects sections defined by sets of beginning and end keywords
long xfdetect_sect(XFILE* f, const kwdset &kwdb, const kwdset &kwde);

/// function resets internal pointer of the given file to the beginning of the actual section
void xf_resetsec(XFILE *f);

/// set section given by enumstr
long xf_setsec(XFILE *f, enumstr sec_alias);

/// set section given by internal section index of the file
long xf_setsec(XFILE *f, long sec_id);

/// copies content of the actual section to the text file
void xf_copysec(XFILE *f, FILE *out);

/// extended fscanf function
int xfscanf(XFILE *in, const char *fmt, ...);

/// returns status of the conversion modifier in the format string
long get_modif(const char *fmt, const char *pconv, const char *modif, int modifl);

/// processes new line characters and line numbers in the xfprintf function
long proc_new_ln_fmt(XFILE *out, const char *fmt);

/// removes comments from the string
void clearcomments(char *lnstr);

/// returns number of conversion in the format string
long checkfmt(const char *fmt);

/// function checks string for integer number
long checkint(const char *lnstr);

/// function checks string for unsigned integer number
long checkuint(const char *lnstr);

/// function checks string for unsigned integer number in octal notation
long checkouint(const char *lnstr);

/// function checks string for unsigned integer number in hexadecimal notation
long checkxuint(const char *lnstr);

/// function checks string for real number
long checkreal(const char *lnstr);

/// function searches assignment suppression character in format string
long check_asterisk(const char *fmt);

// function checks size of the actual line
long check_maxlnsize(XFILE *in, long fbr);

/// function searches conversions in format string
void gettypes(const char *fmt, char *types, const char** ptypes, long nvar);

/// function assembles single conversion from the format string
void getsubfmt(const char *fmt, char *subfmt,long bc);

/// function assembles single conversion from the format string
void getsubfmt(const char *fmt, char *subfmt, char **pconv, long bc);

/// function returns number of located conversions in the format string
long checkfmt(const char *fmt);

/// function removes commented part of string
void clearcomments(char *lnstr);

/// function checks string if it contains only the characters specified in filter - primary for whitespaces 
long isemptystr(const char *lnstr,const char *filter);

/// function detects end of file and eventually prints the error message
void checkfeof(XFILE *in, const char *fmt, long idt, const char *ptype, char *subfmt, char *types, void *apar,
               const char *errfile, int errline, const char *errfunc);

/// checks file for end of section
void checkeos(XFILE *in, const char *fmt, long idt, const char *ptype,char *subfmt, char *types, void *apar, 
              const char *errfile, int errline, const char *errfunc);

/// function checks end of file and end of section
long check_feof_eos(XFILE *in);

/// function searches for errors which occured during xfscanf function and eventually prints the error message
long checkscanferr(XFILE *in, long ret, long reqret, long war, const char *partmsg,
                   const char *fmt, char *subfmt, long idt, const char *ptype, const char *errfile,
                   int errline, const char *errfunc);

void cut_str_sec(XFILE *in, char *lnstr);

/// functions searches for keyword in a string in xfscanf function
long getkwd(XFILE *in,char *lnstr, char *&aptr,char *apar, int &br);

/// functions handles the enum reading (%m conversion) in xfscanf function
long getenum(const char *lnstr, kwdset *akwdset, int *apar, int &br, long &war, int handlingr, long ignorecase, long optenum);

/// function searches for enum keywords in given string 
long matchekwd(const char *ekwd, kwdset *akwdset, int *apar, long ignorecase);

/// function searches for enum integer id in given enum definition akwdset 
long matchekwdint(int *apar, kwdset *akwdset);

/// function searches for composed enum integer id in enum definition akwdset
long matchoptekwdint(int *apar, kwdset *akwdset);

/// function searches for errors which occured during %m conversion in xfscanf function and eventually prints the error message
long checkenumerr(XFILE *in, long ret, long reqret, long war, kwdset *akwdset,
                  const char *fmt, char *subfmt, long idt, const char *ptype, const char *errfile,
                  int errline, const char *errfunc);

/// function finds the first occurence of substring needle in the string haystack
const char *strstr(const char *haystack, const char *needle, long ignorecase);
char *strstr(char *haystack, const char *needle, long ignorecase);

/// function finds the first occurence of substring needle in the string haystack, the case of characters is ignored
char *strstrcis(char *haystack, const char *needle);
const char *strstrcis(const char *haystack, const char *needle);

/// function skips rest of actual line and moves to the next line
void skipline(XFILE *in);

///  function reads the ordinary format characters contained in the format substring fmt from the file in
const char *proc_ord_fmt(XFILE *in, const char *fmt, const char *fmt_end, char *lnstr, char *&astr);

///  function skips the blank characters in the file in
long skip_space(XFILE *in, char *lnstr, char *&astr);

/// function checks whether the '\n' character is included in the accepted set of characters for format "%[...]"
long check_newline_fmt(const char *fmt);

/// function  detects connecting character '@' for multiline string at the end of the of the astr
long detect_multln(char *astr);

/// function determines whether the %a format skips the initial whitespace characters including \n and \r
long a_fmt_skips_whitespaces(const char *fmt);

/// function decomposes the file name to strings with path, filename and suffix
void filename_decomposition (const char *file,char *&path,char *&name,char *&suffix);

/// function decomposes the file name into pointers to the beginnings of path, name and extension
void filename_decomp_ptr (const char *fnstr, const char *&path, const char *&name, const char *&ext);

/// function decomposes the file name into pointers to the beginnings of path, name and extension
void filename_decomp_ptr (char *fnstr, char *&path, char *&name, char *&ext);

/// function returns string with file path from the given filename
char* get_path (char *file, char *&p);
#endif
