#include "stacktrace.h"
#include "systrace.h"

#ifdef TRACE_SOURCE_FILES
 #ifdef LIN_STACK_TRACE
  #include "libtrace.h"
  static char Init_bfd = 1;
 #endif
#endif


/**
  The function initializes internal variable with the program name. It is used
  by Linux in case of usage of bfd library - i.e. macro TRACE_SOURCE_FILES is defined.

  @param fname - string with program name

  @retval The function does not return anything

  Created by Tomas Koudelka, 01.2010
*/
void set_prgname(const char *fname)
{
 #if defined(LIN_STACK_TRACE) && defined(TRACE_SOURCE_FILES)
  set_program_name(fname);
 #endif
}


     
/**
  The function prints stack trace (i.e. call tree) to the text file defined by out.

  @param out - pointer to the opened file
  @param level - level id from which the stack trace will be printed to out

  @return The function does not return anything.

  Created by T. Koudelka, 2008
*/
void stack_trace(FILE *out, long level) 
{
 #ifdef LIN_STACK_TRACE
  lin_stack_trace(out, level);
 #endif
 #ifdef WIN_STACK_TRACE
  win_stack_trace(out, level);
 #endif
}



#ifdef LIN_STACK_TRACE
 #ifndef _GNU_SOURCE
  # define _GNU_SOURCE
 #endif

 #include <stdlib.h>
 #include <dlfcn.h>
 #include <execinfo.h>
 #include <string.h>
 #include <ucontext.h>
 #include <link.h>
 #ifndef NO_CPP_DEMANGLE
  #include <cxxabi.h>
 #endif

 #if defined(REG_RIP)
  # define STACK_IA64
 #elif defined(REG_EIP)
  # define STACK_X86
 #else
  # define STACK_GENERIC
 #endif

/**
  The function prints stack trace (i.e. call tree) to the text file defined by out on Linux system.

  @param out - pointer to the opened file
  @param level - level id from which the stack trace will be printed to out

  @return The function does not return anything.

  Created by T. Koudelka, 2008
*/
void lin_stack_trace(FILE *out, long level) 
{
  ucontext_t *ucontext = new ucontext_t;
  getcontext(ucontext);
  long count = 0;
 #ifdef TRACE_SOURCE_FILES
  char fname[256];
  const char *pname = NULL;
  unsigned int ln;
 #endif

 #if defined(STACK_X86) || defined(STACK_IA64)
  Dl_info dlinfo;
  link_map *plm = NULL;  
  void **bp = 0;
  void *ip = 0;
 #else
  long i;
  void *bt[20];
  char **strings;
  size_t sz;
 #endif

 #if defined(STACK_X86) || defined(STACK_IA64)
  # if defined(STACK_IA64) // 64 bit system
   ip = (void *)ucontext->uc_mcontext.gregs[REG_RIP];
   bp = (void **)ucontext->uc_mcontext.gregs[REG_RBP];
  # elif defined(STACK_X86) // 32 bit system
   ip = (void *)ucontext->uc_mcontext.gregs[REG_EIP];
   bp = (void **)ucontext->uc_mcontext.gregs[REG_EBP];
  # endif

  #ifdef TRACE_SOURCE_FILES
   if ((get_bfd_pointer() == NULL) && (Init_bfd == 1))
   {
     if(dladdr(ip, &dlinfo))
       pname = dlinfo.dli_fname;
     else
       pname = get_program_name();
     if (pname)
     {
       libtrace_init(pname, NULL, NULL);
       Init_bfd = 0;
     }
   } 
  #endif
  while(bp && ip) 
  {
    if(!dladdr1(ip, &dlinfo, (void**)&plm, RTLD_DL_LINKMAP))
      break;

    const char *symname = dlinfo.dli_sname;
   #ifndef NO_CPP_DEMANGLE
    int status;
    char *tmp = abi::__cxa_demangle(symname, NULL, NULL, &status);
    if(status == 0 && tmp)
      symname = tmp;
   #endif

    if (count >= level)
    {
      if (count==level)
      {
        fprintf(out, "\nComplete call trace:\n");
        #ifdef TRACE_SOURCE_FILES
         if (get_bfd_pointer())
         {
           if (plm->l_addr)
             // call for position independent executable => link addres l_addr is not null pointer, dli_fbase is equal to l_addr and
             // instruction pointer ip must be passed relative to dli_fbase address
             libtrace_resolve((void *)((unsigned long)ip - (unsigned long)dlinfo.dli_fbase -1), NULL, 0, fname, 256, ln);
           else
             // call for no_pie gcc option => position dependable executable, dlinfo.dli_fbase is not equal to link address plm->l_addr and
             // instruction pointer is passed directly
             libtrace_resolve((void *)((unsigned long)ip - 1), NULL, 0, fname, 256, ln);

           fprintf(out, "%s [%s, line %u]\n", symname, fname, ln); 
         }
         else
           fprintf(out, "%s (address %p, raddress %#lx)\n", symname, ip, 
                   (unsigned long)(ip) - (unsigned long)(dlinfo.dli_saddr));
        #else
         fprintf(out, "%s (address %p, raddress %#lx)\n", symname, ip, 
                (unsigned long)(ip) - (unsigned long)(dlinfo.dli_saddr));
        #endif
      }
      else
      {
        #ifdef TRACE_SOURCE_FILES
         if (get_bfd_pointer())
         {
           if (plm->l_addr)
             // call for position independent executable => link addres l_addr is not null pointer, dli_fbase is equal to l_addr and
             // instruction pointer ip must be passed relative to dli_fbase address
             libtrace_resolve((void *)((unsigned long)ip - (unsigned long)dlinfo.dli_fbase -1), NULL, 0, fname, 256, ln);
           else
             // call for no_pie gcc option => position dependable executable, dlinfo.dli_fbase is not equal to link address plm->l_addr and
             // instruction pointer is passed directly
             libtrace_resolve((void *)((unsigned long)ip - 1), NULL, 0, fname, 256, ln);

           fprintf(out, "  called by %s [%s, line %u]\n", symname, fname, ln); 
         }
         else
           fprintf(out, "  called by %s (address %p raddress %#lx)\n", symname, ip, 
                   (unsigned long)(ip) - (unsigned long)(dlinfo.dli_saddr));
        #else
         fprintf(out, "  called by %s (address %p raddress %#lx)\n", symname, ip, 
                 (unsigned long)(ip) - (unsigned long)(dlinfo.dli_saddr));
        #endif
      }
    }
   #ifndef NO_CPP_DEMANGLE
    if(tmp)
      free(tmp);
   #endif

    if(dlinfo.dli_sname && !strcmp(dlinfo.dli_sname, "main"))
    {
      fprintf(out, "  in program %s\n\n", dlinfo.dli_fname);
      break;
    }
    count++;
    ip = (void *)bp[1];
    bp = (void **)bp[0];
  }
 
 #else
  sz = backtrace(bt, 20);
  strings = backtrace_symbols(bt, sz);

  fprintf(out, "\nCall trace (for max. %ld levels):\n", sz-level);
  fprintf(out, "%s\n", strings[level]);
  for(i = level+1; i < sz; ++i)
    fprintf(out, "  called by %s()\n", strings[i]);
  fprintf(out, "\n");
 #endif
  delete ucontext;
}
#endif



#ifdef  WIN_STACK_TRACE
 #include <windows.h>
 #ifdef __WATCOMC__
  #include <imagehlp.h>
//  #include "demangle.h"
 #endif
 #ifndef __WATCOMC__
  #include <dbghelp.h>
 #endif

 #define MAX_INSTRUCTION_LENGTH 100
 #define MAX_FUNCTION_PROLOG 100
 #define MAX_SYMBOL_LENGTH 512

// SymSetOptions()
 typedef DWORD (__stdcall *tSSO)(IN DWORD SymOptions);
 tSSO pSSO;

 // SymInitialize()
 typedef BOOL (__stdcall *tSI)(IN HANDLE hProcess, IN PSTR UserSearchPath, IN BOOL fInvadeProcess);
 tSI pSI;

 // SymLoadModule()
 //  typedef DWORD (__stdcall *tSLM)(IN HANDLE hProcess, IN HANDLE hFile, IN PSTR ImageName, 
 //                                      IN PSTR ModuleName, IN DWORD BaseOfDll, IN DWORD SizeOfDll);
 //  tSLM pSLM;
 


#ifdef _X86_
 // SymGetSymFromAddr()
 typedef BOOL (__stdcall *tSGSFA)(IN HANDLE hProcess, IN DWORD dwAddr, OUT PDWORD pdwDisplacement, 
                                  OUT PIMAGEHLP_SYMBOL Symbol);
 tSGSFA pSGSFA;

 // SymGetLineFromAddr()
 typedef BOOL (__stdcall *tSGLFA)(IN HANDLE hProcess, IN DWORD dwAddr,
                                  OUT PDWORD pdwDisplacement, OUT PIMAGEHLP_LINE Line);
 tSGLFA pSGLFA;
#endif

#if defined(_IA64_) || defined(_AMD64_)
 // SymGetSymFromAddr64()
 typedef BOOL (__stdcall *tSGSFA)(IN HANDLE hProcess, IN DWORD64 dwAddr, OUT PDWORD64 pdwDisplacement, 
                                  OUT PIMAGEHLP_SYMBOL64 Symbol);
 tSGSFA pSGSFA;

 // SymGetLineFromAddr64()
 typedef BOOL (__stdcall *tSGLFA)(IN HANDLE hProcess, IN DWORD64 dwAddr,
                                  OUT PDWORD pdwDisplacement, OUT PIMAGEHLP_LINE64 Line);
 tSGLFA pSGLFA;

 // StackWalk64()
 typedef BOOL (__stdcall *tSW)(DWORD MachineType, HANDLE hProcess, HANDLE hThread, LPSTACKFRAME64 StackFrame, 
                               PVOID ContextRecord, PREAD_PROCESS_MEMORY_ROUTINE64 ReadMemoryRoutine,
                               PFUNCTION_TABLE_ACCESS_ROUTINE64 FunctionTableAccessRoutine,
                               PGET_MODULE_BASE_ROUTINE64 GetModuleBaseRoutine,
                               PTRANSLATE_ADDRESS_ROUTINE64 TranslateAddress);
 tSW pSW;

 // SymFunctionTableAccess64()
 typedef PVOID (__stdcall *tSFTA)(HANDLE hProcess, DWORD64 AddrBase);
 tSFTA pSFTA;

 // SymGetModuleBase64()
 typedef DWORD64 (__stdcall *tSGMB)(IN HANDLE hProcess, IN DWORD64 dwAddr);
 tSGMB pSGMB;

#endif

 // UnDecorateSymbolName
 //  typedef DWORD (__stdcall *tUDSN) (IN PCTSTR DecoratedName, OUT PTSTR UnDecoratedName, IN DWORD UndecoratedLength,
 //                                    IN DWORD Flags);
 //  tUDSN pUDSN;

 HINSTANCE hDbghelpDll = NULL;

 
 /**
   The function performs simplified demangling of the function name for WATCOM compiler.
   @param mname - string with mangled name
   @param dname - string for output of demangled name 
                  (the size allocated memory for the output buffer has to be sufficient for the whole of mname string)

   @returns The function writes demangled name in the parameter dname.
   Created by T. Koudelka, 11.2010
 */
 void demangle_watcom(char *mname, char *dname)
 {
   // Standard Watcom prefix name is "W?"
   char *aux = strstr(mname, "W?");
   if (aux) 
   { // prefix was found, we assume the name of the function follows after the "W?"
     aux+=2;
     size_t i = 0;
     // the name of function should be terminated by '$'
     while (aux[i] && aux[i] != '$')
     {
        // copy characters until the end of string or '$' character
        dname[i] = aux[i];
        i++;
     }
     // add string terminating character
     dname[i] = '\0';
   }
   else
     // no mangling prefix was found -> copy the string including terminating \0
     memcpy(dname, mname, strlen(mname)+1);
 }


 /**
   The function returns pointer to the given level of the path specified in the 
   parameter path.

   @param path - string with complete path and file name
   @param level - number of path level which will be left in the returned file name string.
                  Level 1 means that file name only will be copied and no path, level 2 means file name and
                  the name of folder containing the given file, etc. 
                  

   @return The function returns pointer to the file name with reduced path level.

   Created by T. Koudelka, 11.2010
 */
 char *remove_path_level(char *path, long level)
 {
   char *fname = NULL;
   size_t i = strlen(path)-1; // start from the index of last path character
   char c;
   long stop = 0;
   while (path+i >= path)
   {
     c = path[i];
     // backslash or slash path delimiter is checked
     if (c == '\\' || c == '/')
     {
       fname = path+i+1;
       stop++;
       if (stop == level)
         break;
     }
     i--;
   }
   if (fname == NULL)
     fname = path;
   return fname;
 }



/**
  The function prints stack trace (i.e. call tree) to the text file defined by out on Windows system.
  It was tested with MS Visual C++ 2008 Express Edition, Open Watcom C++ 1.9 and Borland C++ Builder 2009.
   * Visual C++ - the switch /clr cannot be used (Common Language Runtime support has to be switched off).
                  No additional steps are necessary. MS native debug format PDB is used by default and
                  it is processed by the dbgHelp function without problems.
   * Watcom C++ - debug format CodeView 4 with ONLY LINE DEBUG SYMBOLS has to be used and stored 
                  in the exe file (not separate symbol file). The settings has to be performed both in
                  target C++ compiler switches and target linker switches. Then the rebase.exe from the 
                  old Windows SDK has to be used in order to strip the debug information out of the exe 
                  file and store them into .dbg file. Use rebase.exe -b 0x400000 -x . prog_name.exe. The 
                  address 0x400000 is the memory segment where exe files are stored by default.
   * C++Builder - debug format CodeView 4 has to be used and stored in the tds file (default storage)
                  The tds2dbg.exe from Bjarne Juul Pasgaard has to be used in order to convert the debug 
                  information from the tds file and store them into .dbg file.

  @param out - pointer to the opened file
  @param level - level id from which the stack trace will be printed to out

  @return The function does not return anything.

  Created by T. Koudelka, 11.2010
*/
 void win_stack_trace(FILE *out, long level)
 {  
   // we load dbghelp.dll dynamically
   if (hDbghelpDll == NULL)
   {
     hDbghelpDll = LoadLibrary(TEXT("dbghelp.dll"));
     if (hDbghelpDll == NULL)
     {
       fprintf(stderr, "\n\nError: Cannot load dbghelp.dll,\n stack trace cannot be performed\n");
       fprintf(stderr, "in file %s, line %d, function %s\n", __FILE__, __LINE__, __FUNCTION__);
       return;
     }

     pSI    = (tSI)    GetProcAddress(hDbghelpDll, "SymInitialize");
     pSSO   = (tSSO)   GetProcAddress(hDbghelpDll, "SymSetOptions");
   #ifdef _X86_
//     pSLM   = (tSLM)   GetProcAddress(hDbghelpDll, "SymLoadModule" );
//     pSW    = (tSW)    GetProcAddress(hDbghelpDll, "StackWalk");
//     pSFTA  = (tSFTA)  GetProcAddress(hDbghelpDll, "SymFunctionTableAccess");
//     pSGMB  = (tSGMB)  GetProcAddress(hDbghelpDll, "SymGetModuleBase");
     pSGSFA = (tSGSFA) GetProcAddress(hDbghelpDll, "SymGetSymFromAddr");
     pSGLFA = (tSGLFA) GetProcAddress(hDbghelpDll, "SymGetLineFromAddr");
//     pUDSN  = (tUDSN)  GetProcAddress(hDbghelpDll, "UnDecorateSymbolName");
   #endif
   #if defined(_IA64_) || defined(_AMD64_)
     pSGSFA = (tSGSFA) GetProcAddress(hDbghelpDll, "SymGetSymFromAddr64");
     pSGLFA = (tSGLFA) GetProcAddress(hDbghelpDll, "SymGetLineFromAddr64");
   #endif

     if ((pSI == NULL)  || (pSSO == NULL) || (pSGSFA == NULL) || (pSGLFA == NULL))
     {
       fprintf(stderr, "Error: some required function not found in imagehlp.dll\n" );
       fprintf(stderr, "in file %s, line %d, function %s\n", __FILE__, __LINE__, __FUNCTION__);
       FreeLibrary(hDbghelpDll);
       return;
     }
   }

   // Initialize symbols
   HANDLE  hprocess = ::GetCurrentProcess ();
   // set options for loading of symbols
   pSSO(SYMOPT_LOAD_LINES| SYMOPT_UNDNAME);
   if ( !pSI(hprocess, NULL, TRUE))
   {
     fprintf(stderr, "Error: debugging symbols cannot be initialized\n");
     fprintf(stderr, "in file %s, line %d, function %s\n", __FILE__, __LINE__, __FUNCTION__);
   }

   // Get current context of processor registers
   CONTEXT context = {0};
   //HANDLE  hthread = ::GetCurrentThread ();
   context.ContextFlags = CONTEXT_FULL;
   #ifdef _X86_
    DWORD reip;
    DWORD rebp;
    DWORD resp;
   #endif

   #if defined(_IA64_) || defined(_AMD64_)
    DWORD64 reip;
    //DWORD64 rebp;
    //DWORD64 resp;
   #endif
   // There is no reliable way for obtaining of the actual context for the 
   // current thread which is portable on the different versions of Windows ->
   // we have to use assembler.
   do{
      #ifdef _X86_
//       __asm {
//         mov rebp, ebp
//         mov resp, esp
//       };
//       context.Ebp = rebp;
//       context.Esp = resp;
       RtlCaptureContext(&context);
      #endif
      #if defined(_IA64_) || defined(_AMD64_)
//        __asm{
//          mov rebp, rbp
//          mov resp, rsp
//        };
//       context.Rbp = rebp;
//       context.Rsp = resp;
       RtlCaptureContext(&context);
      #endif
     // Instruction pointer cannot be captured reliably but 
     // we can use the address of function increased so that the pointer
     // refers somewhere into the code of function defined by 
     // source codes and not to code of function prolog generated by compiler
     // The length of the function prolog is defined by MAX_FUNCTION_PROLOG
     // and it is set to 100.
     #ifdef _X86_
      context.Eip = (unsigned int)(&win_stack_trace) + MAX_FUNCTION_PROLOG;
     #endif
     #if defined(_IA64_) || defined(_AMD64_)
      context.Rip = (long long unsigned int)(&win_stack_trace) + MAX_FUNCTION_PROLOG;
     #endif
   }while(0);

   // Initialize the STACKFRAME structure for the first call. This is only
   // necessary for Intel CPUs, and isn't mentioned in the documentation.
   STACKFRAME sf = { {0} };
   sf.AddrPC.Mode      = AddrModeFlat;
   sf.AddrStack.Mode   = AddrModeFlat;
   sf.AddrFrame.Mode   = AddrModeFlat;

   #ifdef _X86_
    sf.AddrPC.Offset    = context.Eip;
    sf.AddrStack.Offset = context.Esp;
    sf.AddrFrame.Offset = context.Ebp;
   #endif
   #if defined(_IA64_) || defined(_AMD64_)
    sf.AddrPC.Offset    = context.Rip;
    sf.AddrStack.Offset = context.Rsp;
    sf.AddrFrame.Offset = context.Rbp;
   #endif

   // Walk the stack up to DATA_LEVELS calls
   int count=0;
   reip = 0;
   IMAGEHLP_LINE Line;
   memset(&Line, 0, sizeof(Line));
   Line.SizeOfStruct = sizeof(Line);
  #if defined(_IA64_) || defined(_AMD64_)
   // Get a stack trace for this frame by StackWalk64 function i.e. pSW (works under Visual Studio only)
   while (::pSW(IMAGE_FILE_MACHINE_I386, GetCurrentProcess(), GetCurrentThread(), 
          &sf, &context, NULL, pSFTA, pSGMB, NULL))
  #else
   while (sf.AddrFrame.Offset != 0)
  #endif
   {
     // IMAGEHLP is wacky, and requires you to pass in a pointer to an
     // IMAGEHLP_SYMBOL structure. The problem is that this structure is
     // variable length. That is, you determine how big the structure is
     // at runtime. This means that you can't use sizeof(struct).
     // So...make a buffer that's big enough, and make a pointer
     // to the buffer. We also need to initialize not one, but TWO
     // members of the structure before it can be used.
     BYTE symbol_buffer[sizeof(IMAGEHLP_SYMBOL) + MAX_SYMBOL_LENGTH];
     PIMAGEHLP_SYMBOL psymbol = (PIMAGEHLP_SYMBOL)symbol_buffer;
     memset (symbol_buffer, 0, sizeof (symbol_buffer));
     psymbol->SizeOfStruct = sizeof(symbol_buffer);
     psymbol->MaxNameLength = 512;

     // structure used for line number and source file name retrieving 
     IMAGEHLP_LINE Line;
     memset(&Line, 0, sizeof(Line));
     Line.SizeOfStruct = sizeof(Line);
     // auxiliary pointer to the source file name with reduced length of path
     char *fname = NULL;
     DWORD offset_l = 0;
     #ifdef _X86_
      DWORD sym_displacement = 0;// Displacement of the input address, relative to the start of the symbol
     #endif
     #if defined(_IA64_) || defined(_AMD64_)
      DWORD64 sym_displacement = 0;// Displacement of the input address, relative to the start of the symbol
     #endif

     // the first condition is used in the case, that the loop is controlled
     // by the StackWalk function and the second is for direct registry acces which has 
     // to be used for Watcom compiler
     if ((sf.AddrPC.Offset == sf.AddrReturn.Offset) || (sf.AddrPC.Offset == reip))
     {
       fprintf(out, "\nError: Endless loop in the stack trace\n");
       break;
     }
         
     // backup of return point of the actual stack level
     reip = sf.AddrPC.Offset;

     // Get symbols from the actual address by function SymGetSymFromAddr i.e. pSGSFA
     if (pSGSFA != NULL)
     {
       if (pSGSFA(::GetCurrentProcess(), sf.AddrPC.Offset, &sym_displacement, psymbol) == FALSE)
         fprintf(out, "\nCannot get the symbols: %u\n", (unsigned int)(GetLastError()));
     }
     // show line number info, NT5.0-method (SymGetLineFromAddr) i.e. pSGLFA
     if (pSGLFA != NULL)
     { // SymGetLineFromAddr was loaded
         
       // Watcom supports only CodeView format version 4 which has to be stored 
       // in the separated file progname.dbg, where progname is the name of the executed program.
       // The .DBG file can be created by rebase.exe utility from the Windows SDK (1999) which
       // extracts the debug information stored in the exe file
       // In the case of CV4 format, the line numbers and source file names can be detected  
       // only for the instruction address of the line beginning and not for arbitrary addresses  
       // between the beginning of the call instruction and the return point.
       // Thus we have to test addresses starting one byte before return point until the
       // line beginning of the function call is captured.
       // The Visual Studio uses PDB format which accepts arbitrary instruction address in the range 
       // belonging to the call instruction of the given function. In such the case, the loop finishes
       // after the first pass.
       unsigned range = 1; // starting range from the address of return point
       while ((pSGLFA(::GetCurrentProcess(), sf.AddrPC.Offset-range, &offset_l, &Line) == FALSE) && (range < MAX_INSTRUCTION_LENGTH))
       {
         // increase range
         range++;
       }
       if (pSGLFA(::GetCurrentProcess(), sf.AddrPC.Offset-range, &offset_l, &Line) != FALSE)
       { // line and filename were located
         // remove whole path except of the first folder preceding the source file
         fname = remove_path_level(Line.FileName, 2);
       }
     }
         // buffer for the demangled name of the function
     char tname[MAX_SYMBOL_LENGTH];
     tname[0] = '\0';
     // start log from the required stack frame level 
     if (count >= level)
     {
       // demangling of the function name
       #ifdef __WATCOMC__
        // demangling function from the Watcom library - prepared for the future use (proper header files have to be added)
//       __demangle_l(psymbol->Name, strlen(psymbol->Name), tname, MAX_SYMBOL_LENGTH);
        // simplified demangling function
        demangle_watcom(psymbol->Name, tname);
       #else
//        DWORD err;
//        // demangler for Visual Studio is not necessary, the symbol names are demangled by
//        // SYMOPT_UNDNAME option during initiaton phase (pSI)
//        err = pUDSN(psymbol->Name, tname, strlen(psymbol->Name)+1, UNDNAME_NAME_ONLY);
        memcpy(tname, psymbol->Name, strlen(psymbol->Name)+1);
       #endif
       if (count==level)
       {
         fprintf(out, "\nComplete call trace:\n");
         fprintf(out, "%s ", tname);
         if (fname == NULL)
           fprintf(out, "(address %p, raddress %#lx)\n", (void*)sf.AddrPC.Offset, 
                   (unsigned long)(sf.AddrPC.Offset) - (unsigned long)(psymbol->Address));
         else
           fprintf(out, "[%s, line %lu]\n", fname, Line.LineNumber);
       }
       else
       {
         fprintf(out, "  called by %s ", tname);
         if (fname == NULL)
           fprintf(out, "(address %p, raddress %#lx)\n", (void*)sf.AddrPC.Offset, 
                   (unsigned long)(sf.AddrPC.Offset) - (unsigned long)(psymbol->Address));
         else
           fprintf(out, "[%s, line %lu]\n", fname, Line.LineNumber);
       }
     }
     if (strcmp(psymbol->Name, "main") == 0)
       break;
     count++;
     // Unwind the stack one level up
     #ifdef _X86_
      // return point is stored before ebp address 
      sf.AddrPC.Offset    = *((DWORD *)(sf.AddrFrame.Offset+sizeof(void*))); 
      // the content of ebp address is base pointer for the upper stack level
      sf.AddrFrame.Offset = *((DWORD *)(sf.AddrFrame.Offset));        
      // stack pointer address esp can be set to arbitrary value above the base pointer ebp
      // (it is not required but the Windows function could require a valid address)
      sf.AddrStack.Offset = sf.AddrFrame.Offset-sizeof(void*);
     #endif
     #if defined(_IA64_) || defined(_AMD64_)
      // return point is stored before ebp address 
      // sf.AddrPC.Offset    = *((DWORD64 *)(sf.AddrFrame.Offset+sizeof(void*))); 
      // the content of ebp address is base pointer for the upper stack level
      // sf.AddrFrame.Offset = *((DWORD64 *)(sf.AddrFrame.Offset));        
      // stack pointer address esp can be set to arbitrary value above the base pointer ebp
      // (it is not required but the Windows function could require a valid address)
      // sf.AddrStack.Offset = sf.AddrFrame.Offset-sizeof(void*);
     #endif
   }
   fprintf(out, "\n");
 }
#endif
