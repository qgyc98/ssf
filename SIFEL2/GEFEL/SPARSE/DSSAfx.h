// DSSAfx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//
#if !defined(_DSSAFX_H__)
#define _DSSAFX_H__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

// support for Watcom 10.6
#ifdef __WATCOMC__
 #if __WATCOMC__ < 1100
  #ifndef BOOL_DEF
   enum bool {false = 0, true = 1};
   #define BOOL_DEF
  #endif 
 #endif
#endif

// This is just for compilation with MFC
#ifdef BUILDING_CCFEModel
   #include <afx.h>
   #ifdef _DEBUG
      #define new DEBUG_NEW
      #undef THIS_FILE
      static char THIS_FILE[] = __FILE__;
   #endif
#endif

//#include <omp.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>


#ifndef min2
#define min2(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef max2
#define max2(a,b) ((a)>(b)?(a):(b))
#endif


typedef int                 BOOL;
typedef unsigned char       BYTE;
typedef unsigned long ULONG;

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef PURE 
#define PURE =0
#endif//PURE 

#ifndef UNREFERENCED_PARAMETER
#define UNREFERENCED_PARAMETER(P) (P)
#endif//UNREFERENCED_PARAMETER

#define EXIT_NOT_IMPLEMENTED printf("Not implemented method called! - exiting"); exit(-1);


#define DSS_NAMESPASE_BEGIN    
#define DSS_NAMESPASE_END 
//#define DSS_NAMESPASE_BEGIN namespace DSS {
//#define DSS_NAMESPASE_END }

#endif // !defined(_DSSAFX_H__)

