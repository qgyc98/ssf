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

#ifndef _TREE
#define _TREE

#include <stddef.h>
#include "func.h"

union NumVar
{
  int WhichVar;
  double Number;
};

struct Node
{
  Node *Parent;
  Node *Child_One;
  Node *Child_Two;
  func_type Type;
  NumVar Value;
  int Brackets;
  double CurrentValue;
  double (*evalfn)(double);  
};

class Parser;

class TTree
{
 public:
  Node *Root;

  TTree();
  ~TTree();
  void InitNode(Node *, func_type);
  void InitNode(Node *, func_type, NumVar, int brack=0);
  Node *AddParent(Node *);
  Node *AddChildOne(Node *);
  Node *AddChildTwo(Node *);
  Node *DeleteNode(Node *,Node *Blank=NULL,Node *Other=NULL);
  void Reset(int SetTo,Node *From=NULL);
  Node *Copy(Node *From,Node *To);
  void Delete(Node *From);
  Node *AddBelowOne(Node *, func_type, double Value=0);
  Node *AddBelowTwo(Node *, func_type, double Value=0);

};

#endif
