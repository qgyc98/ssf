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

/* Header file for Parser class */

#ifndef _TPARSER
#define _TPARSER

#include "equ.h"
#include "func.h"

struct Token
{
   Token *Last;
   Token *Next;
   func_type Type;
   NumVar Value;
};

class Parser
{
public:
  Parser();
  ~Parser();
  Equation *TextToTree(char *);
  void TreeToText(Equation *);
  void TidyUpEqn(Equation *);

private:
  void PassOne();
  void PassTwo();
  void ErrorCheck();
  void Digit(char *, double *, bool *, int *);
  void Alpha(const char *,char *);
  void Symbol(char *);
  void EndFunc(char *);
  void EndNum(double *, bool *, int *);
  void AddToken(func_type);
  bool Priority(func_type, func_type);
  void RemoveBrackets();
  void TidyUp();
  void TidyNode(Node *&);
  void TreeToText();
  char * TreeToText(Node *);
  void NodeText(int, Node *, char *, int &);

private:
  char *EquationText;
  Token *EquationTokenStart;
  Token *EquationTokenEnd;
  Equation *equation;
};

#endif
