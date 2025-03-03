/*
 * mtk - Maths Toolkit for X11
 *
 * Copyright 1994-1997   andrewr@chiark.greenend.org.uk (Andrew Ross)
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

#ifndef _MOREMATH
#define _MOREMATH

#define EINEXACT  	100
#define EDIVZERO 	101

double factorial(double);
double combination(double,double);
double permutation(double,double);
double integ(double);
double frac(double);
double step(double);
double sec(double);
double cosec(double);
double cot(double);

double ansi_asinh(double);
double ansi_acosh(double);
double ansi_atanh(double);
double ansi_lgamma(double);
double ansi_jn(int,double);
double ansi_erf(double);
double ansi_erfc(double);

#endif
