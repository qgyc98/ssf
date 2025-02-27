/* 
 * Originally, derived by Aurelian Melinte from addr2line.c and associated binutils files, version 2.18. 
 * rewritten by Tomas Koudelka
 */


#include "systrace.h"

#ifdef TRACE_SOURCE_FILES
#ifdef BFD_VER
 #if BFD_VER <= 233
  #define BFD_GET_SECTION_FLAGS(abfd, asection)  bfd_get_section_flags((abfd), (asection))
  #define BFD_GET_SECTION_VMA(abfd, asection)    bfd_get_section_vma((abfd), (asection))
  #define BFD_GET_SECTION_SIZE  bfd_get_section_size
 #else
  #define BFD_GET_SECTION_FLAGS(abfd, asection)  bfd_section_flags((asection))
  #define BFD_GET_SECTION_VMA(abfd, asection)    bfd_section_vma((asection))
  #define BFD_GET_SECTION_SIZE  bfd_section_size 
 #endif
#endif

#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include "iotools.h"

//#pragma GCC visibility push(hidden) 
#include <bfd.h>
#define DMGL_PARAMS	 (1 << 0)	/* Include function args */
#define DMGL_ANSI	 (1 << 1)	/* Include const, volatile, etc */
//#pragma GCC visibility pop 



typedef struct _libtrace_data {
    bfd_boolean unwind_inlines;  /* Unwind inlined functions. */
    bfd_boolean with_functions;  /* Show function names.  */
    bfd_boolean do_demangle;     /* Demangle names.  */
    bfd_boolean base_names;      /* Strip directory names.  */
    asymbol **syms;              /* Symbol table.  */

    bfd *abfd;                   /* Binary File Descriptor*/
    asection *section;           /* */
} libtrace_data;



static libtrace_data m_libtrace_data = {
    FALSE, // unwind_inlines
    FALSE, // with_functions
    FALSE, // do_demangle
    TRUE,  // base_names
    NULL,  // symbol table
    NULL,  // bfd
    NULL,  // section
};



/** 
  These variables are used to pass information between
  translate_addresses and find_address_in_section.  
*/
typedef struct _sym_info {
    bfd_vma pc;
    const char *filename;
    const char *functionname;
    unsigned int line;
    bfd_boolean found;
} sym_info; 



static int slurp_symtab (bfd *);
static void find_address_in_section (bfd *, asection *, void *);
static void find_offset_in_section (bfd *, asection *, sym_info *);
static int translate_addresses (bfd *abfd, asection *section, 
                                void *addr, 
                                char *buf_func, size_t buf_func_len, 
                                char *buf_file, size_t buf_file_len,
                                unsigned int &line);

/// name of the executable file
static const char *program_name = NULL; 



/**
  The function sets internal variable program_name.

  @param name - pointer to string with the program name
  
  @retval The function does not return anything.

  Created by Tomas Koudelka, 01.2010
*/
void set_program_name(const char *name)
{
  program_name = name;
}



/**
  The function returns pointer to string with the program name.
  The string is stored in the internal variable program_name.

  @retval pointer to string with the program name

  Created by Tomas Koudelka, 01.2010
*/
const char *get_program_name()
{
  return program_name;
}



void *get_bfd_pointer()
{
  return (void *)m_libtrace_data.abfd;
}



/**
  The function prints an error caused during bfd initialization.
  
  @param string - string with user defined part of error message
  
  @retval The function does not return anything.
*/
void bfd_nonfatal(const char *string)
{
  const char *errmsg = bfd_errmsg (bfd_get_error ());
  
  if (string)
    print_err("%s: %s", __FILE__, __LINE__, __func__, string, errmsg);
  else
    print_err("%s", __FILE__, __LINE__, __func__, errmsg);
}



/**
  The function returns the size of the named file.  If the file does not
  exist, or if it is not a real file, then a suitable error message is printed 
  and zero is returned.
  
  @param file_name - string with file name

  @retval 0 - in case of failure
  @retval file size in bytes in case of success 
   
  Created by Tomas Koudelka, 01.2010
*/
off_t get_file_size(const char * file_name)
{
  struct stat statbuf;
    
  if (stat(file_name, &statbuf) < 0) 
  {
    if (errno == ENOENT)
      print_err("'%s': no such file", __FILE__, __LINE__, __func__, file_name);
    else
      print_warning("could not locate '%s'.  reason: %s",
                    __FILE__, __LINE__, __func__, file_name, strerror (errno));
  } 
  else 
  {
    if (! S_ISREG (statbuf.st_mode))
      print_warning("'%s' is not an ordinary file", __FILE__, __LINE__, __func__, file_name);
    else
      return statbuf.st_size;
  }
  
  return 0;
}



/**
  The function prints the possible matching targets after a FALSE return 
  from bfd_check_format_matches with bfd_get_error () == bfd_error_file_ambiguously_recognized.

  @param p - pointer to array of strings with target names terminated by NULL pointer

  @retval The function does not return anything.

  Created by Tomas Koudelka, 01.2010
*/
void list_matching_formats(char **p)
{
  if (!p || !*p)
    return;

  char **x = p;
  char *buf;
  char *aux;
  long l = 0;

  while (*p)
  {
    l += strlen(*p)+1;
    p++;
  }
  p = x;
  buf = new char[l];
  aux = buf;
  while (*p)
  {
    sprintf(aux, " %s", *p);
    aux += strlen(*p)+1;
    p++;
  }
  print_err("matching formats:%s", __FILE__, __LINE__, __func__, buf);
  delete [] buf;
}



/**
   The function sets the default BFD target based on the configured target.  Doing
   this permits the binutils to be configured for a particular target,
   and linked against a shared BFD library which was configured for a
   different target.*/
/*
void set_default_bfd_target(void)
{
  // The macro TARGET is defined by Makefile. 
  //   E.g.: -DTARGET='"i686-pc-linux-gnu"'.
  const char *target = TARGET;
  
  if (! bfd_set_default_target (target)) 
  {
    print_err("can't set BFD default target to `%s': %s",
              __FILE__, __LINE__, __func__, target, bfd_errmsg (bfd_get_error ()));
    return;
  }
  return; 
}
*/



/** 
  The function reads the symbol table.

  @param abfd - pointer to basic data structure of binary file descriptor 

  @retval 0 - on success
  @retval -1 - in case of an failure

  Created by Tomas Koudelka, 01.2010
*/
static int slurp_symtab(bfd *abfd)
{
  long symcount;
  unsigned int size;
  
  if ((bfd_get_file_flags(abfd) & HAS_SYMS) == 0)
    return -1;
  
  symcount = bfd_read_minisymbols(abfd, FALSE, (void **)(&m_libtrace_data.syms), &size);
  if (symcount == 0)
    symcount = bfd_read_minisymbols(abfd, TRUE /* dynamic */, 
                                    (void **)(&m_libtrace_data.syms), &size);
  
  if (symcount < 0) 
  {
    bfd_nonfatal(bfd_get_filename(abfd));
    return -1;
  }
    
  return 0; 
}



/**
  The function looks for an code address in a section.  This is called via
  bfd_map_over_sections.  

  @param abfd - pointer to basic data structure of binary file descriptor
  @param section - pointer to the required section
  @param data - pointer to user data, in this case it requires structure which holds source file name, 
                line number and function name
*/
static void find_address_in_section(bfd *abfd, asection *section, void *data)
{
  bfd_vma vma;
  bfd_size_type size;
  sym_info *psi = (sym_info*)data; 
  
  if (psi->found)
    return;
  
  if ((BFD_GET_SECTION_FLAGS(abfd, section) & SEC_ALLOC) == 0)
    return;
  
  vma = BFD_GET_SECTION_VMA(abfd, section);
  if (psi->pc < vma)
    return;
  
  size = BFD_GET_SECTION_SIZE(section);
  if (psi->pc >= vma + size)
    return;
  
  psi->found = bfd_find_nearest_line(abfd, section, 
                                     m_libtrace_data.syms, psi->pc - vma,
                                     &psi->filename, &psi->functionname, 
                                     &psi->line);
}



/** 
  The function looks for an offset in a section. This is directly called.

  @param abfd - pointer to basic data structure of binary file descriptor
  @param section - pointer to the required section
  @param psi - pointer to structure which holds source file name, line number and function name

  @retval The function does not return anything.

  Created by Tomas Koudelka, 01.2010
*/
static void find_offset_in_section(bfd *abfd, asection *section, sym_info *psi)
{
  bfd_size_type size;
  
  if (psi->found)
    return;
  
  if ((BFD_GET_SECTION_FLAGS(abfd, section) & SEC_ALLOC) == 0)
    return;
  
  size = BFD_GET_SECTION_SIZE(section);
  if (psi->pc >= size)
    return;
  
  psi->found = bfd_find_nearest_line(abfd, section, 
                                     m_libtrace_data.syms, psi->pc,
                                     &psi->filename, &psi->functionname, 
                                     &psi->line);
}



/**
  The function translates the given code address to the source file name and line (optionally also function name
  if set with_functions flag in structure libtrace_data).
  
  @param bfd - pointer to basic data structure of binary file descriptor
  @param section - pointer to the required section
  @param addr - void pointer with required hexadecimal address code
  @param buf_func - buffer for the function name (can be NULL)
  @param buf_func_len - length of the buf_func, no more characters then this value is written to this buffer
  @param buf_file - buffer for the source file name
  @param buf_file_len - length of the buf_file, no more characters then this value is written to this buffer
  @param line - address line number in the buf_file

  @retval 0 - on success
  
  Created by Tomas Koudelka, 01.2010
*/
static int translate_addresses(bfd *abfd, asection *section, 
                               void *xaddr, 
                               char *buf_func, size_t buf_func_len, 
                               char *buf_file, size_t buf_file_len,
                               unsigned int &line)  
{
  #define ADDR_BUF_LEN ((CHAR_BIT/4)*(sizeof(void*))+1)
  char addr[ADDR_BUF_LEN+1] = {0};
  sym_info si = {0, NULL, NULL, 0, false}; 
  
  sprintf(addr, "%p", xaddr);
  si.pc = bfd_scan_vma (addr, NULL, 16);

  si.found = FALSE;
  if (section)
    find_offset_in_section(abfd, section, &si);
  else
    bfd_map_over_sections(abfd, find_address_in_section, &si);

  if (! si.found) 
  {
    if (buf_func != NULL)
      snprintf(buf_func, buf_func_len, "%s ??:0", 
               m_libtrace_data.with_functions ? "??" : "");
    if (buf_file != NULL)
      snprintf(buf_file, buf_file_len, "%s ??:0", 
               m_libtrace_data.with_functions ? "??" : "");
    line = 0;
  } 
  else 
  {
    line = si.line;
    do 
    {
      if (m_libtrace_data.with_functions) 
      {
        const char *name;
        char *alloc = NULL;
          
        name = si.functionname;
        if (name == NULL || *name == '\0')
          name = "??";
        else 
          if (m_libtrace_data.do_demangle) 
          {
            alloc = bfd_demangle(abfd, name, DMGL_ANSI | DMGL_PARAMS);
            if (alloc != NULL)
              name = alloc;
          }

        if (buf_func != NULL)
          snprintf(buf_func, buf_func_len, "%s", name);

        if (alloc != NULL)
          free (alloc);
      }

      if (m_libtrace_data.base_names && si.filename != NULL) 
      {
        const char *h = strrchr(si.filename, '/');
        if (h != NULL)
          si.filename = h + 1;
      }

      if (buf_file != NULL)
        snprintf(buf_file, buf_file_len, "%s", 
                 si.filename ? si.filename : "??");
      if (!m_libtrace_data.unwind_inlines)
        si.found = FALSE;
      else
        si.found = bfd_find_inliner_info(abfd, &si.filename, &si.functionname, &si.line);
    } while (si.found);
  }
  
  return si.found;
}



/**
  The function initializes system of detection of source file name and line
  for given program addresses.

  @param file_name    - name of the executable program
  @param section_name - name of searched section (NULL can be used as default value)
  @param target       - format name of the executable program (NULL can be used as default value)

  @retval  0 - on success
  @retval -1 - in case of failure

  Created by Tomas Koudelka 01.2010
*/
int libtrace_init(const char *file_name, const char *section_name, const char *target)
{
  char **matching = NULL;
  
  bfd_init();
  //set_default_bfd_target();
  
  if (get_file_size(file_name) < 1)
    return -1;
  
  m_libtrace_data.abfd = bfd_openr(file_name, target);
  if (m_libtrace_data.abfd == NULL) 
  {
    bfd_nonfatal(file_name);
    return -1; 
  }
  
  if (bfd_check_format(m_libtrace_data.abfd, bfd_archive)) 
  {
    print_err ("%s: cannot get addresses from archive",
               __FILE__, __LINE__, __func__, file_name);
    return -1; 
  }
  
  if (!bfd_check_format_matches(m_libtrace_data.abfd, bfd_object, &matching)) 
  {
    bfd_nonfatal(bfd_get_filename(m_libtrace_data.abfd));
    if (bfd_get_error() == bfd_error_file_ambiguously_recognized) 
    {
      list_matching_formats(matching);
      free (matching);
    }
    return -1;
  }
  
  if (section_name != NULL) 
  {
    m_libtrace_data.section = bfd_get_section_by_name(m_libtrace_data.abfd, 
                                                      section_name);
    if (m_libtrace_data.section == NULL)
      print_err("%s: cannot find section %s",
                __FILE__, __LINE__, __func__, file_name, section_name);
  } else
    m_libtrace_data.section = NULL;
  
  if (0 != slurp_symtab(m_libtrace_data.abfd))
    return -1;
  
  return 0; 
}



/**
  The function closes system of detection of source file name and line
  for given program addresses and it frees memory allocated for this system.

  @retval The function does not return anything.

  Created by Tomas Koudelka 01.2010
*/
void libtrace_close(void)
{
  if (m_libtrace_data.syms != NULL) 
  {
    free (m_libtrace_data.syms);
    m_libtrace_data.syms = NULL;
  }
  
  bfd_close(m_libtrace_data.abfd);
}



/**
  The function returns source file name and line (optionally also function name if set 
  with_functions flag in structure libtrace_data) for the given code address.
  
  @param addr - void pointer with required hexadecimal address code
  @param buf_func - buffer for the function name (can be NULL)
  @param buf_func_len - length of the buf_func, no more characters then this value is written to this buffer
  @param buf_file - buffer for the source file name
  @param buf_file_len - length of the buf_file, no more characters then this value is written to this buffer
  @param line - address line number in the buf_file

  @retval 0 - on success
  
  Created by Tomas Koudelka, 01.2010
*/
int libtrace_resolve(void *addr, 
                     char *buf_func, size_t buf_func_len, 
                     char *buf_file, size_t buf_file_len,  
                     unsigned int &line)
{
  int ret = FALSE; 
  ret = translate_addresses(m_libtrace_data.abfd, 
                            m_libtrace_data.section, addr, 
                            buf_func, buf_func_len,
                            buf_file, buf_file_len, line);
  assert(0 == ret); 
  
  return 0; 
}



/*
  Sample of usage of the above functions


  libtrace_init(file_name, NULL, NULL)) 

  void  *sym = NULL;     // required code address should be obtained somewhere else
  char  file[PATH_MAX+1] = {0}; 
  unsigned int  line;
  .
  .
  .
  sym = ......
  libtrace_resolve(sym, NULL, 0, file, PATH_MAX, line);
  .
  .
  .
  libtrace_close();
*/

#endif

