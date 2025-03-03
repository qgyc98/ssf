#ifndef DEF_VECTOR
#define DEF_VECTOR
#include <stdio.h>
#include <assert.h>
#include <string.h>

//#define DEBUG_VECTOR

#include "iotools.h"
#include "stackdef.h"

struct vecstat
{
  unsigned dealloc:1;  ///< flag for vector memory deallocation (0=no deallocation, 1=release memory by delete)
  vecstat() {dealloc=0;};
};


struct vector
/**
   This file declares struct of the vector, which implements
   %vector with elements type of double.
   There is also declarations of the functions for the %vector computing.
*/
{
  long n;    ///< the number of %vector entries (components)
  double *a; ///< pointer to memory ,where are stored elements of %vector
  long size; ///< real length of array a (due to memory reallocation)
  vecstat stat; ///< bit array of vector status flags (deallocation)

  inline vector() { n = size = 0L; a = NULL;}; ///< default constructor
  vector(long n);               ///< allocating constructor
  vector(const vector &v);      ///< copy constructor

  /**
    The constructor assigns preallocated memory to the member array a.
   
    @param n is the number of %vector elements
    @param dealloc is the flag for deallocation control (0=no deallocation ,1=release memory by delete operator)
    @param ptr is the pointer to the allocated memory which will be used for storage of n elements
   
    created  4.5.2015 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline vector (long n, unsigned dealloc, double *ptr)
  {
   #ifdef STCKLIM  
    if (ptr == NULL)
      print_err("cannot allocate memory on stack for vector v(%ld)", __FILE__, __LINE__, __func__, n);
   #endif

    vector::n = size = n;
    a = ptr;
    stat.dealloc = dealloc;

    memset (a, 0, n*sizeof(*a));

   #ifdef STCKLIM  
    if (dealloc == 0)
    {
      Sma += size;
      if (Sma > Smm)
        Smm = Sma;
    }
   #endif

   #ifdef DEBUG_VECTOR
    Acv++;
    Ava += size;
    if (Ava > Avmax)
      Avmax = Ava;
   #endif
  };


  /**
     The function operator enables access to the elements of the member array a
     
     @param i is the number of the desired %vector element
     @retval Returns i-th element of %vector
     
     created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
     modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline double& operator () (long i)
  {
    #ifdef DEBUG_VECTOR
    if ((i >= n) || (i < 0))
       print_err("vector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };
  
  /**
     The function operator enables CONSTANT access to the elements of the member array a
     
     @param i is the number of the desired %vector element
     @retval Returns i-th element of %vector
     
     created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
     modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline double& operator () (long i) const
  {
    #ifdef DEBUG_VECTOR
    if ((i >= n) || (i < 0))
       print_err("vector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };
  
  /**
     The array index operator enables access to the elements of the member array a
     
     @param i is the number of the desired %vector %element
     @retval i-th %element of %vector
     
     created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
     modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline double& operator [] (long i)
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       print_err("vector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };
  
  /**
     The array index operator enables CONSTANT access to the elements of the member array a
     
     @param i is the number of the desired %vector element
     @retval i-th %element of %vector
     
     created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
     modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline double& operator [] (long i) const
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       print_err("vector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };

  double get_item(long i) { return a[i]; }
  void set_item(long i, double val) { a[i] = val; }
  
  inline ~vector() ///< destructor
  {
   #ifdef DEBUG_VECTOR
    if (a != NULL)
    {
      Acv--;
      Ava -= size;
    }
   #endif

    if (stat.dealloc)
      delete [] a;
   #ifdef STCKLIM  
    else
      Sma -= size;
   #endif
  };
};


struct ivector
/**
   This file declares struct of the vector, which implements
   %vector with %elemnts type of long.
*/
{
  long n;   ///< number of %vector elements
  long *a;  ///< pointer to memory ,where are stored elements of %vector
  long size; ///< real length of array a (due to memory reallocation)
  vecstat stat; ///< bit array of vector status flags (deallocation)
  
  inline ivector() { n = size = 0L; a = NULL;}; ///< default constructor
  ivector(long n);                      ///< allocating constructor
  ivector(const ivector &v);            ///< copy constructor
  ivector& operator=(const ivector &v); ///< copy assignment operator

  /**
    The constructor assigns preallocated memory to the member array a.
   
    @param n is the number of %vector elements
    @param dealloc is the flag for deallocation control (0=no deallocation ,1=release memory by delete operator)
    @param ptr is the pointer to the allocated memory which will be used for storage of n elements
   
    created  4.5.2015 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline ivector (long n, unsigned dealloc, long *ptr)
  {
   #ifdef STCKLIM  
    if (ptr == NULL)
      print_err("cannot allocate memory on stack for vector v(%ld)", __FILE__, __LINE__, __func__, n);
   #endif

    ivector::n = size = n;
    a = ptr;
    stat.dealloc = dealloc;

    memset (a, 0, n*sizeof(*a));
  
   #ifdef STCKLIM  
    if (dealloc == 0)
    {
      Sma += size;
      if (Sma > Smm)
        Smm = Sma;
    }
   #endif

   #ifdef DEBUG_VECTOR
    Acv++;
    Ava += size;
    if (Ava > Avmax)
      Avmax = Ava;
   #endif
  };
  
  /**
     The function operator enables access to the elements of the member array a
     
     @param i is the number of the desired %vector element
     @retval Returns i-th element of %vector
     
     created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
     modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline long& operator () (long i)
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       print_err("ivector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };

  /**
     The function operator enables CONSTANT access to the elements of the member array a
     
     @param i is the number of the desired %vector element
     @retval Returns i-th element of %vector
     
     created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
     modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline long& operator () (long i) const
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       print_err("ivector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };

  /**
   The array index operator enables access to the elements of the member array a
   
   @param i is the number of the desired %vector %element
   @retval i-th %element of %vector
   
   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline long& operator [] (long i)
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       print_err("ivector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };

  /**
   The array index operator enables CONSTANT access to the elements of the member array a
   
   @param i is the number of the desired %vector %element
   @retval i-th %element of %vector
   
   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline long& operator [] (long i) const
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       print_err("ivector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };

  inline ~ivector() ///< destructor
  {
   #ifdef DEBUG_VECTOR
    if (a != NULL)
    {
      Aciv--;
      Aiva -= size;
    }
   #endif

    if (stat.dealloc)
      delete [] a;
   #ifdef STCKLIM  
    else
      Sma -= size;
   #endif
  };


  long get_item(long i) { return a[i]; }
  void set_item(long i, long val) { a[i] = val; }
};

/**
   this structure serves for defined attributes
   it is used e.g. in stochastic analysis
   
   TKo
*/
struct atsel{
  ///  number of components
  long num;
  ///  number of basic components (optional parameter)
  long nba;
  ///  number of internal (hidden) components (optional parameter)
  long nia;

  long *atrib;

  atsel (){
    num=0;  
    nba=0;  
    nia=0;  
    atrib=NULL;
  };
  ~atsel (){  delete [] atrib;  };

  long read(XFILE* in);
};

#ifdef DEBUG_VECTOR
 long give_acv();
 long give_avmax();
#endif

/// allocates %vector
long allocv(long n, vector &vec);
/// allocates %ivector
long allocv(long n, ivector &vec);
/// reallocates %vector
long reallocv(long n, vector &vec);
/// reallocates %vector with help of preallocated memory
long reallocv(long n, vector &vec, unsigned dealloc, double *ptr);
/// creates %vector which contains reference to %vector src
long makerefv(vector &ref, vector &src);
/// creates %vector which contains reference to array ptr
long makerefv(vector &ref, double *ptr, long n);
/// creates %vector which contains reference to %vector src
long makerefv(ivector &ref, ivector &src);
/// creates %vector which contains reference to array ptr
long makerefv(ivector &ref, long *ptr, long n);
/// reallocates %ivector
long reallocv(long n, ivector &vec);
/// reallocates %ivector with help of preallocated memory
long reallocv(long n, ivector &vec, unsigned dealloc, long *ptr);

/// copies contents of %vector
long copyv(const vector &src, vector &dest);
/// copies contents of %ivector
long copyv(const ivector &src, ivector &dest);
/// copies contents of array
long copyv(const long *src,long *dest,long n);
/// copies contents of array
long copyv(const double *src, double *dest, long n);
/// copies contents of array to vector
long copyv(const double *src, vector &dest);
/// copies contents of array to vector
long copyv(const long *src, ivector &dest);
/// copies contents of vector to array
long copyv(const vector &src, double *dest);
/// copies contents of ivector to array
long copyv(const ivector &src, long *dest);

/// swaps pointers to component arrays of two vectors, i.e. swaps their content
long swapv(vector &a, vector &b);

/// copies range of components between vectors
long rcopyv(const vector &src, long src_fi, vector &dest, long dest_fi, long n);
/// copies range of components between vectors given by arrays
long rcopyv(const double *src, long src_fi, long src_n, double *dest, long dest_fi, long dest_n, long n);

/// copies contents of %vector multiplied by scalar
long copymultv(const vector &src, vector &dest, double c);
/// copies contents of %ivector multiplied by scalar
long copymultv(const ivector &src, ivector &dest, long c);
/// copies contents of array multiplied by scalar
long copymultv(const double *src, double *dest, double c, long n);

/// fills contents of %vector with value c
long fillv(double c, vector &vec);
/// fills contents of %ivector with value c
long fillv(long c, ivector &vec);
/// fills contents of integer %vector given by array with value c
void fillv(long c, long *vec,long n);
/// fills contents of %vector given by array with value c
void fillv(double c, double *vec,long n);

/// fills array with zeros
long nullv (double *a,long n);
/// fills integer array with zeros
long nullv (long *a,long n);
/// fills %vector with zeros
long nullv (vector &a);
/// fills %ivector with zeros
long nullv (ivector &a);

/// changes sign of %vector components
void chsgnv(vector &v);
/// changes sign of array components
void chsgnv(double *a, long n);

/// deallocates %vector
long destrv(vector &vec);
/// deallocates %ivector
long destrv(ivector &vec);

/// adds 2 vectors c = a+b
long addv(const vector &a, const vector &b, vector &c);
/// adds vectors a = a+b
long addv(vector &a, const vector &b);
/// adds 2 double arrays, 2nd array is multplied by constant c = a + bc*b
void addmultv(const vector &a, const vector &b, double bc, vector &c);
/// adds 2 double arrays, 2nd array is multplied by constant a = a + bc*b
void addmultv(vector &a, const vector &b, double bc);

/// adds 2 ivectors c = a+b
long addv(const ivector &a, const ivector &b, ivector &c);

/// adds 2 double arrays a = a+b
void addv(double *a,double *b,long n);

/// adds 2 double arrays c = a+b
void addv(double *a,double *b, double *c, long n);

/// adds 2 double arrays, 2nd array is multplied by constant and added to the first array
void addmultv (double *a, double *b, double bc, long n);
/// adds 2 double arrays, 2nd array is multplied by constant, the result is stored in the third array
void addmultv (double *a, double *b, double bc, double *c, long n);

/// adds 2 vectors, both are multiplied by the independent constants and stored in the resulting vector
void addmultv (vector &a, double ac, vector &b, double bc, vector &c);
/// adds 2 double arrays, both are multiplied by the independent constants and stored in the first array
void addmultv (double *a, double ac, double *b, double bc, long n);
/// adds 2 double arrays, both are multiplied by the independent constants and stored in the last array
void addmultv (double *a, double ac, double *b, double bc, double *c, long n);

/// subtracts 2 vectors
long subv(const vector &a, const vector &b, vector &c);
/// subtracts 2 vectors
long subv(vector &a, const vector &b);
/// subtracts 2 ivectors
long subv(const ivector &a, const ivector &b, ivector &c);

/// subtracts 2 arrays
void subv(double *a,double *b,long n);
/// subtracts 2 arrays, result is stored to the another array
void subv (double *a,double *b,double *c,long n);

/// performs cross product of 2 vectors
long crprd(const vector &a, const vector &b, vector &c);
/// performs cross product of 2 ivectors
long crprd(const ivector &a, const ivector &b, ivector &c);

/// performs scalar product of 2 vectors
long scprd(const vector &a, const vector &b, double &scalar);
/// performs scalar product of 2 vectors
double scprd(const vector &a, const vector &b);
/// performs scalar product of 2 ivectors
long scprd(const ivector &a, const ivector &b, long &scalar);
/// performs scalar product of 2 arrays
double scprd(const double *a, const double *b, long n);
/// performs scalar product of 2 vectors
double scprd_t(const vector &a, const vector &b, long nt);
/// performs scalar product of two arrays on several threads
double scprd_t(const double *a, const double *b, long n, long nt);
/// performs scalar product of 2 arrays
double ss (double *a,double *b,long n);

/// multiplies %vector by constant
long cmulv(double c, const vector &u, vector &v);
/// multiplies %ivector by constant of type double
long cmulv(double c, const ivector &u, vector &v);
/// multiplies %ivector by constant of type long
long cmulv(long c, const ivector &u, ivector &v);

/// multiplies %vector by constant
long cmulv(double c, vector &a);
/// multiplies %ivector by constant of type long
long cmulv(long c, ivector &a);
/// multiplies double array by constant
void cmulv(double c, double *a, long n);
/// multiplies double array by constant
void cmulv(double c, double *a, double *u, long n);

/// computes norm of %vector
double normv(const vector &a);
/// computes norm of %vector
double normv(const double *a, long n);
/// computes norm of %ivector
double normv(const ivector &a);

/// computes selective norm of %vector
double snormv(const vector &a, long fi, long nc);
/// computes selective norm of %vector
double snormv(const double *a, long n, long fi, long nc);

/// normalizes %vector components so that after function call |a| = 1
void normalize(vector &a);
/// normalizes %vector components and stores them into b, |b| = 1
void normalize(const vector &a, vector &b);
/// normalizes %vector components so that after function call |a| = 1
void normalize(double *a, long n);
/// normalizes %vector components and stores them into b, |b| = 1
void normalize(double *a, long n, double *b);

/// returns maximum of vector components
//double maxcompv(const vector &a);
//double maxcompv(const double *a, long n);

/// cosine angle of 2 vectors
long cosav(const vector &a, const vector &b, double &cos);
/// cosine angle of 2 ivectors
long cosav(const ivector &a, const ivector &b, double &cos);

/// absolute value of vector components
void absv(vector &a);
/// absolute value of vector components
void absv(ivector &a);

/// prints contents of %vector to the file out with given precision and field width
long printv(const vector &u, FILE *out = stdout, int prec = 3, int width = 11);
/// prints contents of %vector given by array to the file out with given precision and field width
long printv(const double *u, long n, FILE *out = stdout, int prec = 3, int width = 11);
/// prints %vector u to the file out
long printv(FILE *out, vector &u);
/// prints contents of %ivector
long printv(const ivector &u, FILE *out = stdout, int width = 11);
/// prints contents of integer %vector given by array
long printv(const long *u, long n, FILE *out = stdout, int width = 11);

/// reads contents of %vector from file
long readv(XFILE *in, vector &u);
/// reads contents of %vector from string
long readv(char *in, vector &u);
/// reads contents of integer %vector from file
long readv(XFILE *in, ivector &u);
/// reads contents of integer %vector from string
long readv(char *in, ivector &u);

/// function extracts ncomp components (from index fi) from %vector b
void extract (vector &a, const vector &b, long fi, long ncomp);

/// extracts only positive components from %vector a and stores them to %vector b
long extractposv (const vector &a, vector &b);
/// extracts only negative components from %vector a and stores them to %vector b
long extractnegv (const vector &a, vector &b);

/// Shell sort of ivector
void shell_sort(ivector &v);

/// Shell sort for long integer array
void shell_sort(long *v, long n);

#endif
