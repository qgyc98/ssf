/**
  The file defines/undefines macros controlling stacktrace and
  libtrace depending on various operating systems and installed
  system libraries.
 
  __linux__ and __LINUX__ are macros defined by compilers in case of Linux target

  _WIN32 and __WINDOWS__ are macros defined by compilers in case of Windows target

  LIN_STACK_TRACE is user defined macro controlling stack trace under Linux

  WIN_STACK_TRACE is user defined macro controlling stack trace under Windows

  TRACE_SOURCE_FILES is user defined macro controlling output of source file names and lines under Linux.
    If defined, auxiliary bfd  library (binutils-dev) have to be installed on the Linux system.
*/

// Predefined macros in GCC
#ifdef __linux__
 #define LIN_STACK_TRACE
// #define TRACE_SOURCE_FILES
#endif

// Predefined macros in Open Watcom compiler
#ifdef __LINUX__
 #define LIN_STACK_TRACE
 #define TRACE_SOURCE_FILES
#endif

// Predefined macros in Visual C++, Borland C++ and Open Watcom compilers
#ifdef _WIN32
 #define WIN_STACK_TRACE
#endif

// Predefined macros in Watcom compiler
#ifdef __WINDOWS__
 #define WIN_STACK_TRACE
#endif
