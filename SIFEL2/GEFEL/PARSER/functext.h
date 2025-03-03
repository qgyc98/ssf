/*
 * mtk - Maths Toolkit for X11
 *
 * Copyright 1994-1996   andrewr@chiark.greenend.org.uk (Andrew Ross)
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 ********/

/* Header file for function text */

#ifndef _FUNCTEXT
#define _FUNCTEXT

#include <math.h>
#include "moremath.h"
#include "func.h"

typedef struct mathfunc_s
{
  const char *func_text;
  double (*func_ptr)(double);
} mathfunc;

#ifdef PARSER_CPP

mathfunc func_details[MTK_NUMBER_OF_FUNCS+1] = {
  {"",NULL},  /* "" */
  {" + ",NULL},  /* + */
  {" - ",NULL},  /* - */
  {" * ",NULL},  /* * */
  {" / ",NULL},  /* / */
  {"C",NULL},  /* C */
  {"P",NULL},  /* P */
  {"^",NULL},  /* ^ */
  {"ln ",log},
  {"log ",log10},
  {"exp ",exp},
  {"sin ",sin},
  {"cos ",cos},
  {"tan ",tan},
  {"sec ",sec},
  {"cosec ",cosec},
  {"cot ",cot},
  {"arcsin ",asin},
  {"arccos ",acos},
  {"arctan ",atan},
  {"sinh ",sinh},
  {"cosh ",cosh},
  {"tanh ",tanh},
  {"arsinh ",ansi_asinh},
  {"arcosh ",ansi_acosh},
  {"artanh ",ansi_atanh},
  {"erf ",ansi_erf},
  {"erfc ",ansi_erfc},
  {"bessj ",NULL},
  {"H ",step},
  {"abs ",fabs},
  {"int ",integ},
  {"frac ",frac},
  {"!",factorial},
  {"(",NULL},
  {"",NULL},
  {"",NULL},
  {"pi ",NULL},
  {")",NULL},
  {",",NULL}
};

#else

extern mathfunc func_details[39];

#endif

#endif /* _FUNCTEXT */
