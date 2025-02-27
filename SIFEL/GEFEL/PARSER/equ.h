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

/* Header file for Equation class */

#ifndef _EQU
#define _EQU

#define DEGTORAD   0.01745329252
#define EZERO      0

#include "tree.h"
#include "qlist.h"

struct range
{
  double From;
  double To;
  double Step;
  double Accuracy;
};

struct variable
{
  char Name[12];
  double Value;
};

class Parser;

// support for Watcom 10.6 
#ifdef __WATCOMC__
 #if __WATCOMC__ < 1100
  #ifndef BOOL_DEF
   enum bool {false = 0, true = 1};
   #define BOOL_DEF
  #endif
 #endif 
#endif

class Equation : public TTree
{
  friend class Parser;

public:
  Equation();
  ~Equation();
  double Evaluate();
  void Differentiate( Equation *, int, Parser * );
  double NumInt( range, int );
  double FindRoot( range, int );
  double RK4( range, int, int );
  void PowerSeries( Equation *, int, int, Parser * );
  const char *VarText(int);
  double VarValue(int);

private:
  void RemoveBrackets();
  void EvalNode(Node *);
  void DiffFromNode(Node *, Node **);
  variable *NumToVarPtr(int);
  void AddPowerSeriesTerm(int, Node *, double, double, Parser *);

public:
  QList<variable> Variables;
  char *EquationText;
  int DiffVar;
  bool Radians;

private:
  int brackets;
};

#endif
