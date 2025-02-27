#ifndef STACKDEF_H
#define STACKDEF_H

#if defined(__linux__)  || defined(__APPLE__)  // Linux or Apple GCC compiler
 #define DEFALLOCA <alloca.h>
 #define ALLOCAF alloca
#else  // Open Watcom, Borland C++, MS Visual C++
 #define DEFALLOCA <malloc.h>

 #ifdef _MSC_VER  // MS Visual C++
  #define ALLOCAF _alloca
 #else
  #define ALLOCAF alloca // Open Watcom, Borland C++
 #endif
#endif

//
// uncomment the following line to switch back to the original heap allocation version of %vector and %matrix classes
//
//#undef DEFALLOCA


#ifdef DEFALLOCA
 #include DEFALLOCA
 /** 
   Limit number of double elements allocated in stack memory by alloca() 
   function called in the allocations of %vector, %ivector and %matrix class instances.
   The total ammount of memory allowed for the dynamic allocation is given
   by STCKLIM*sizeof(double). 

   If the STCKLIM is defined then the memory allocated is checked in order not to exceed 
   the given limit otherwise the checks are off and the fastest allocation is performed.
 */
// #define STCKLIM 150000  // comment out this line to obtain the fastest stack allocation
#endif


#ifndef STCKGLOBVAR
 // this case is common for every source or header file except of vector.cpp
 extern long Sma; ///< actual number of double elements allocated on stack for all instances of %vector, %ivector and %matrix classes
 extern long Smm; ///< maximum number of double elements allocated on stack in the allocation history of all %vector, %ivector and %matrix instances
#else
 // this case is invoked in the vector.cpp (due to global variable intialization)
 long Sma=0; ///< actual number of double elements allocated on stack for all instances of %vector, %ivector and %matrix classes
 long Smm=0; ///< maximum number of double elements allocated on stack in the allocation history of all %vector, %ivector and %matrix instances
#endif


#ifdef DEFALLOCA // alloca() function is used for dynamic allocation in the stack memory

 #ifdef STCKLIM // all stack memory allocation is checked -> slower but 'safe'
  /** 
    Macro for allocation of %vector in stack memory which is intended
    to use in %vector class one parametric constructors and
    allocv function that will be repleced automatically by their three/four 
    parametric variants by the macro expansion. If the STACKLIM would 
    be exceeded by the actual stack memory allocation then the memory 
    is allocated from the heap.
    The macro call ASTCKVEC(n) is expanded into '(n), dealloc, (ptr)' 
    where (n) is the macro argument representing the number of allocated 
    elements, (ptr) is pointer to the memory allocated either on the 
    stack or heap and dealloc is 0 (stack memory) or 1 (heap memory)
  */  
  #define ASTCKVEC(n)\
          (n),\
          ((Sma + (n) <= STCKLIM) ? 0 : 1),\
          ((Sma + (n) <= STCKLIM) ? (double*)ALLOCAF((n)*sizeof(double)) : new double[(n)])

  /** 
    Macro for reallocation of %vector memory which is intended
    to use in connection with two parametric reallocv function 
    that will be repleced automatically by its four parametric variant 
    by the macro expansion. Reallocation of either stack or heap memory
    is supported with respect to stack memory usage. If the STCKLIM would 
    exceeded by the actual stack memory allocation then the memory is 
    allocated from heap.
    The macro call RSTCKVEC(n) is expanded into '(n), (v), (dealloc), (ptr)' 
    where (n) is the macro argument representing the number of allocated 
    elements, (v) is the reference to the %vector reallocated, (ptr) is pointer 
    to the memory allocated either on the stack or heap if the memory reallocation 
    is required or NULL if no reallocation is necessary and dealloc is 0 
    (reallocation of stack memory/ nor reallocation) or 1 (reallocation of heap memory)
  */
  #define RSTCKVEC(n , v)\
          (n),\
          (v),\
          (((n) > (v).size) ? ((Sma + (n) <= STCKLIM) ? 0 : 1) : 0),\
          (((n) > (v).size) ? ((Sma + (n) <= STCKLIM) ? (double*)ALLOCAF((n)*sizeof(double)) : new double[(n)]) : NULL)

  /** 
    Macro for allocation of %ivector in stack memory which is intended
    to use in %ivector class single parametric constructors and
    allocv function that will be repleced automatically by their three/four 
    parametric variants by the macro expansion. If the STACKLIM would 
    be exceeded by the actual stack memory allocation then the memory 
    is allocated from the heap.
    The macro call ASTCKIVEC(n) is expanded into '(n), dealloc, (ptr)' 
    where (n) is the macro argument representing the number of allocated 
    elements, (ptr) is pointer to the memory allocated either on the 
    stack or heap and dealloc is 0 (stack memory) or 1 (heap memory)
  */  
  #define ASTCKIVEC(n)\
          (n),\
          ((Sma + (n) <= STCKLIM) ? 0 : 1),\
          ((Sma + (n) <= STCKLIM) ? (long*)ALLOCAF((n)*sizeof(long)) : new long[(n)])

  /** 
    Macro for reallocation of %ivector memory which is intended
    to use in connection with two parametric reallocv function 
    that will be repleced automatically by its four parametric variant 
    by the macro expansion. Reallocation of either stack or heap memory
    is supported with respect to stack memory usage. If the STCKLIM would 
    exceeded by the actual stack memory allocation then the memory is 
    allocated from heap.
    The macro call RSTCKIVEC(n) is expanded into '(n), (v), dealloc, (ptr)' 
    where (n) is the macro argument representing the number of allocated 
    elements, (v) is the reference to the %ivector reallocated, (ptr) is 
    the pointer to the memory allocated either on the stack or heap if 
    the memory reallocation is required or NULL if no reallocation is 
    necessary and dealloc is 0 (reallocation of stack memory/ nor reallocation) 
    or 1 (reallocation of heap memory)
  */
  #define RSTCKIVEC(n , v)\
          (n),\
          (v),\
          (((n) > (v).size) ? ((Sma + (n) <= STCKLIM) ? 0 : 1) : 0),\
          (((n) > (v).size) ? ((Sma + (n) <= STCKLIM) ? (long*)ALLOCAF((n)*sizeof(long)) : new long[(n)]) : NULL)

  /** 
    Macro for allocation of %matrix in stack memory which is intended
    to use in %matrix/%imatrix class two parametric constructors and allocm function 
    that will be repleced automatically by their four/five parametric 
    variants by the macro expansion. If the STCKLIM would be exceeded by 
    the actual stack memory allocation then the memory is allocated from 
    the heap.
    The macro call ASTCKMAT(m, n) is expanded into '(m), (n), dealloc, (ptr)' 
    where (m) and (n) are the macro argument representing the number of allocated 
    %matrix rows and columns respectively, (ptr) is pointer to the memory allocated 
    either on the stack or heap and dealloc is 0 (stack memory) or 1 (heap memory)
  */  
  #define ASTCKMAT(m, n)\
          (m),\
          (n),\
          ((Sma + ((m)*(n)) <= STCKLIM) ? 0 : 1),\
          ((Sma + ((m)*(n)) <= STCKLIM) ? (double*)ALLOCAF((m)*(n)*sizeof(double)) : new double[(m)*(n)])

  /** 
    Macro for reallocation of %matrix memory which is intended
    to use in connection with three parametric reallocm function 
    that will be repleced automatically by its five parametric variant 
    by the macro expansion. Reallocation of either stack or heap memory
    is supported with respect to stack memory usage. If the STCKLIM would 
    exceeded by the actual stack memory allocation then the memory is 
    allocated from heap.
    The macro call RSTCKMAT(n) is expanded into '(m), (n), (b), dealloc, (ptr)' 
    where (m) and (n) are the macro argument representing the number of allocated 
    %matrix rows and columns respectively, (b) is the reference to the %matrix reallocated, 
    (ptr) is the pointer to the memory allocated either on the stack or heap if the memory 
    reallocation is required or NULL if no reallocation is necessary and dealloc 
    is 0 (reallocation of stack memory/ nor reallocation) or 1 (reallocation of heap memory)
  */
  #define RSTCKMAT(m, n , b)\
          (m),\
          (n),\
          (b),\
          ((((m)*(n)) > (b).size) ? ((Sma + ((m)*(n)) <= STCKLIM) ? 0 : 1) : 0),\
          ((((m)*(n)) > (b).size) ? ((Sma + ((m)*(n)) <= STCKLIM) ? (double*)ALLOCAF((m)*(n)*sizeof(double)) : new double[(m)*(n)]) : NULL)

  /** 
    Macro for allocation of identity %matrix in stack memory which is intended
    to use in allocim function that will be repleced automatically by its four parametric 
    variant by the macro expansion. If the STCKLIM would be exceeded by the actual stack 
    memory allocation then the memory is allocated from the heap.
    The macro call ASTCKIMAT(n, b) is expanded into '(n), (b), dealloc, (ptr)' 
    where (m) and (n) are the macro argument representing the number of allocated 
    %matrix rows and columns respectively, (b) is the %matrix reallocated, (ptr) is 
    the pointer to the memory allocated either on the stack or heap and dealloc 
    is 0 (stack memory) or 1 (heap memory)
  */  
  #define ASTCKIMAT(n, b)\
          (n),\
          (b),\
          ((Sma + ((n)*(n)) <= STCKLIM) ? 0 : 1),\
          ((Sma + ((n)*(n)) <= STCKLIM) ? (long*)ALLOCAF((n)*(n)*sizeof(long)) : new long[(n)*(n)])









 #else  // no stack limit is given -> the fastest allocation but more dangerous

  /** 
    Macro for allocation of %vector in stack memory which is intended to use 
    as argument wrapper in %vector class one parametric constructors 
    and allocv function that will be repleced automatically by their three/four 
    parametric variants by the macro expansion.
    The macro call ASTCKVEC(n) is expanded into '(n), 0, (ptr)' 
    where (n) is the macro argument representing the number of allocated 
    elements, (ptr) is pointer to the memory allocated on the stack.
  */  
  #define ASTCKVEC(n)\
          (n),\
          0,\
          ((double*)ALLOCAF((n)*sizeof(double)))

  /** 
    Macro for reallocation of %vector memory from stack which is intended
    to use as argument wrapper of  two parametric reallocv function that will 
    be repleced automatically by its four parametric variant by the 
    macro expansion. 
    The macro call RSTCKVEC(n, v) is expanded into '(n), (v), 0, (ptr)' 
    where (n) is the macro argument representing the number of allocated 
    elements, (v) is the reference to the given reallocated %vector, (ptr) 
    is pointer to the memory allocated on the stack if the reallocation is 
    required or NULL if no reallocation is necessary.
  */
  #define RSTCKVEC(n , v)\
          (n),\
          (v),\
          0,\
          (((n) > (v).size) ? (double*)ALLOCAF((n)*sizeof(double)) : NULL)

  /** 
    Macro for allocation of %ivector in stack memory which is intended to use 
    as argument wrapper in %ivector class single parametric constructors 
    and allocv function that will be repleced automatically by their three/four 
    parametric variants by the macro expansion.
    The macro call ASTCKIVEC(n) is expanded into '(n), 0, (ptr)' 
    where (n) is the macro argument representing the number of allocated 
    elements, (ptr) is pointer to the memory allocated on the stack.
  */  
  #define ASTCKIVEC(n)\
          (n),\
          0,\
          ((long*)ALLOCAF((n)*sizeof(long)))

  /** 
    Macro for reallocation of %ivector memory from stack which is intended
    to use as argument wrapper of  two parametric reallocv function that will 
    be repleced automatically by its four parametric variant by the 
    macro expansion. 
    The macro call RSTCKIVEC(n, v) is expanded into '(n), (v), (ptr), 0' 
    where (n) is the macro argument representing the number of allocated 
    elements, (v) is the reference to the given reallocated %ivector, (ptr) 
    is pointer to the memory allocated on the stack if the reallocation is 
    required or NULL if no reallocation is necessary.
  */
  #define RSTCKIVEC(n , v)\
          (n),\
          (v),\
          0,\
          (((n) > (v).size) ? (long*)ALLOCAF((n)*sizeof(long)) : NULL)

  /** 
    Macro for allocation of %matrix in stack memory which is intended
    to use in %matrix/%imatrix class two parametric constructors and allocm function 
    that will be repleced automatically by their four/five parametric 
    variants by the macro expansion.
    The macro call ASTCKMAT(m, n) is expanded into '(m), (n), 0, (ptr)' 
    where (m) and (n) are the macro argument representing the number of allocated 
    %matrix rows and columns respectively and (ptr) is pointer to the memory allocated 
    on the stack.
  */  
  #define ASTCKMAT(m, n)\
          (m),\
          (n),\
          0,\
          ((double*)ALLOCAF((m)*(n)*sizeof(double)))

  /** 
    Macro for reallocation of %matrix memory which is intended
    to use in connection with two parametric reallocm function 
    that will be repleced automatically by its five parametric variant 
    by the macro expansion. Reallocation of either stack or heap memory
    is supported with respect to stack memory usage. If the STCKLIM would 
    exceeded by the actual stack memory allocation then the memory is 
    allocated from heap.
    The macro call RSTCKMAT(n) is expanded into '(m), (n), (b), (ptr), 0' 
    where (m) and (n) are the macro argument representing the number of allocated 
    %matrix rows and columns respectively, (b) is the %matrix reallocated, (ptr) is 
    the pointer to the memory allocated on the stack if the memory reallocation 
    is required or NULL if no reallocation is necessary.
  */
  #define RSTCKMAT(m, n , b)\
          (m),\
          (n),\
          (b),\
          0,\
          ((((m)*(n)) > (b).size) ? (double*)ALLOCAF((m)*(n)*sizeof(double)) : NULL)

  /** 
    Macro for allocation of identity %matrix in stack memory which is intended
    to use in allocim function that will be repleced automatically by its four parametric 
    variant by the macro expansion.
    The macro call ASTCKIMAT(n, b) is expanded into '(n), (b), (ptr), 0' 
    where (n) is the macro argument representing the number of allocated 
    %matrix rows and columns, (b) is the reference to the %matrix reallocated and 
    (ptr) is pointer to the memory allocated on the stack.
  */  
  #define ASTCKIMAT(n, b)\
          (n),\
          (b),\
          0,\
          ((long*)ALLOCAF((n)*(n)*sizeof(long)))

 #endif









#else  // alloca() function is not available -> safe but slow heap allocations will be used

 /** 
   Macro for allocation of %vector in stack memory which is intended
   to use in %vector class one parametric constructors and
   allocv function that will be called with heap memory allocations.  
   The macro call ASTCKVEC(n) is simply expanded into '(n)' where 
   (n) is the number of %vector components.
 */  
  #define ASTCKVEC(n) (n)

  /** 
    Macro for reallocation of %vector memory  which is intended
    to use as argument wrapper of  %vector class two parametric 
    reallocv function that will be called with heap memory allocations.
    The macro call RSTCKVEC(n, v) is simply expanded into '(n), (v)' 
    where (n) is the number of %vector components and (v) is the reference
    to the %vector reallocated.
  */
  #define RSTCKVEC(n, v) (n), (v)

 /** 
   Macro for allocation of %ivector in stack memory which is intended
   to use in %ivector class one parametric constructors and
   allocv function that will be called with heap memory allocations.  
   The macro call ASTCKIVEC(n) is simply expanded into '(n)' where 
   (n) is the number of %ivector components.
 */  
  #define ASTCKIVEC(n) (n)

  /** 
    Macro for reallocation of %ivector memory  which is intended
    to use as argument wrapper of  %ivector class two parametric 
    reallocv function that will be called with heap memory allocations.
    The macro call RSTCKIVEC(n, v) is simply expanded into '(n), (v)'
    where (n) is the number of %ivector components and (v) is the reference
    to the %ivector reallocated.  
  */
  #define RSTCKIVEC(n, v) (n), (v)

  /** 
    Macro for allocation of %matrix/%imatrix in stack memory which is intended
    to use in %matrix class one parametric constructor and
    allocm function that will be called with heap memory allocations.  
    The macro call ASTCKMAT(m, n) is simply expanded into '(m), (n)'
    where (m) is the number of rows and (n) is the number of columns.
  */  
  #define ASTCKMAT(m, n) (m), (n)

  /** 
    Macro for reallocation of %matrix memory  which is intended
    to use as argument wrapper of two parametric reallocm function 
    that will be called with heap memory allocations.
    The macro call RSTCKMAT(m, n, b) is simply expanded into '(m), (n), (b)'
    where (m) is the number of rows, (n) is the number of columns and
    (b) is the reference to the %matrix reallocated.
  */
  #define RSTCKMAT(m, n, b)   (m), (n), (b)

  /** 
    Macro for allocation of identity %matrix memory which is intended
    to use as argument wrapper of two parametric allocim function 
    that will be called with heap memory allocations.
    The macro call ASTCKIMAT(n, b) is simply expanded into '(n), (b)'
    where (n) is the number of rows or columns and (b) is the reference 
    to the %matrix allocated.
  */
  #define ASTCKIMAT(n, b)   (n), (b)
 #endif


#endif
