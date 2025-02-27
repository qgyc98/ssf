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

/* Header file for function and error constants */

#ifndef _FUNC
#define _FUNC

typedef int func_type;

#define MTK_NOFUNC      0
#define MTK_ADD		1
#define MTK_SUBTRACT	2
#define MTK_MULTIPLY	3
#define MTK_DIVIDE	4
#define MTK_COMBINATION	5
#define MTK_PERMUTATION	6
#define MTK_POWER	7
#define MTK_LN		8
#define MTK_LOG		9
#define MTK_EXP		10
#define	MTK_SIN		11
#define	MTK_COS		12
#define	MTK_TAN		13
#define MTK_SEC		14
#define MTK_COSEC	15
#define MTK_COT		16
#define MTK_ARCSIN	17
#define MTK_ARCCOS	18
#define MTK_ARCTAN	19
#define MTK_SINH	20
#define MTK_COSH	21
#define MTK_TANH	22
#define MTK_ARSINH	23
#define MTK_ARCOSH	24
#define MTK_ARTANH	25
#define MTK_ERF         26
#define MTK_ERFC        27
#define MTK_BESSJ       28
#define MTK_STEP        29
#define MTK_ABS		30
#define MTK_INT		31
#define MTK_FRAC	32
#define MTK_FACTORIAL	33
#define MTK_LBRACKET	34
#define MTK_NUMBER	35
#define MTK_VARIABLE	36
#define MTK_PI_SYM	37
#define MTK_RBRACKET	38
#define MTK_COMMA       39
#define MTK_NUMBER_OF_FUNCS 39


#endif /* _FUNC */

