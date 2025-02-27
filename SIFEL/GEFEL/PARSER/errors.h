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

/* Header file for error types */

#ifndef _ERRORS
#define _ERRORS

enum error_type {
NO_ERROR,
ERR_OPERATOR,
ERR_BRACKET,
ERR_NO_TEXT,
ERR_FUNCTION,
ERR_DOM,
ERR_RANGE,
ERR_VARIABLES,
ERR_DIV_ZERO,
ERR_INEXACT,
ERR_DIFF,
ERR_NO_ROOT,
ERR_GRAPH,
ERR_SYMBOL,
NUM_ERRS
};

#ifdef ERROR_MSGS

const char *error_msgs[NUM_ERRS] = {
  "No error",                      // NO_ERROR
  "Illegal use of operator",       // ERR_OPERATOR
  "Incorrect use of brackets",     // ERR_BRACKET
  "No equation entered",           // ERR_NO_TEXT
  "Illegal use of function",       // ERR_FUNCTION
  "Function argument out of range",  // ERR_DOM
  "Function value out of range",   // ERR_RANGE
  "Undefined variable",            // ERR_VARIABLES
  "Divide by zero error",          // ERR_DIV_ZERO
  "Loss of precision",             // ERR_INEXACT
  "Unable to differentiate function",   // ERR_DIFF
  "Are you sure there is a root in the interval",   // ERR_NO_ROOT
  "Unable to graph function",      // ERR_GRAPH
  "Unknown symbol in equation"     // ERR_SYMBOL
};

#endif /* ERROR_MSGS */

#endif /* _ERRORS */
