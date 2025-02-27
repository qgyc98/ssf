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

// Equation class

#include <math.h>
//#include <bits/nan.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include "moremath.h"
#include "equ.h"
#include "parser.h"
#include "func.h"
#include "functext.h"
#include "errors.h"

//Constructor for the Equation class
Equation::Equation() : TTree()
{
  Radians = true;
  EquationText = new char[1];
  EquationText[0] = 0;
  Variables.setAutoDelete(1);
}

//Destructor for the Equation class
Equation::~Equation()
{
  delete [] EquationText;
  Variables.clear();
}

// Evaluate function
double 
Equation::Evaluate()
{
  Node *TreePos;
  variable *Variable;
  TreePos=Root;
  errno = 0;//set system variable to zero to prevent crashes in parser from other preceding routines
	
  Reset(0);
  try {
    while (Root->Brackets==0) {
      while (TreePos->Child_One!=NULL) {
        if ((TreePos->Child_One)->Brackets==1) break;
        TreePos=TreePos->Child_One;
      };
      switch (TreePos->Type) {
      case MTK_NUMBER :
        (TreePos->CurrentValue)=(TreePos->Value).Number;
        (TreePos->Brackets)=1;
        break;
      case MTK_VARIABLE :
        Variable=Variables.at( TreePos->Value.WhichVar - 1);
        if (Variable==NULL) {
          throw ERR_VARIABLES;
        };
        TreePos->CurrentValue = Variable->Value;
        TreePos->Brackets = 1;
        break;
      case MTK_PI_SYM :
        TreePos->CurrentValue = 3.14159265358979323846;
        TreePos->Brackets = 1;
        break;
      default :
        break;
      };
      TreePos = TreePos->Parent;
      if (TreePos==NULL) break;
      if (TreePos->Child_Two==NULL) {
        EvalNode(TreePos);
      }
      else {
        if (TreePos->Child_Two->Brackets==0) {
          TreePos=TreePos->Child_Two;
        }
        else {
          EvalNode(TreePos);
        };
      };
    };
  }
  catch(...) {
    Reset(2);
    throw;
  }
  return (Root->CurrentValue);
}

// Evaluate from the given node
void 
Equation::EvalNode(Node *TreePos)
{
  if (errno) {
    switch (errno) {
    case EDIVZERO :
      errno=EZERO;
      throw ERR_DIV_ZERO;
    case ERANGE :
      errno=EZERO;
      throw ERANGE;
    case EINEXACT :
      errno=EZERO;
      throw ERR_INEXACT;
    };
  };
 
  if ( TreePos->evalfn ) {
    if ( Radians || (TreePos->Type < MTK_SIN) || ( TreePos->Type > MTK_ARCTAN) ) {
      TreePos->CurrentValue = (*(TreePos->evalfn))((TreePos->Child_One)->CurrentValue);
    }
    else {
      if ( TreePos->Type <= MTK_COT ) 
        TreePos->CurrentValue = (*(TreePos->evalfn))((TreePos->Child_One)->CurrentValue*DEGTORAD);
      else
        TreePos->CurrentValue = (*(TreePos->evalfn))((TreePos->Child_One)->CurrentValue)/DEGTORAD;
    }
  }
  else {
    switch (TreePos->Type) {
    case MTK_ADD :
      TreePos->CurrentValue=((TreePos->Child_One)->CurrentValue)
        + ((TreePos->Child_Two)->CurrentValue);
      break;
    case MTK_SUBTRACT :
      TreePos->CurrentValue=((TreePos->Child_One)->CurrentValue)
        - ((TreePos->Child_Two)->CurrentValue);
      break;
    case MTK_MULTIPLY :
      TreePos->CurrentValue=((TreePos->Child_One)->CurrentValue)
        * ((TreePos->Child_Two)->CurrentValue);
      break;
    case MTK_DIVIDE :
      if (TreePos->Child_Two->CurrentValue == 0) {
        throw ERR_DIV_ZERO;
      };
      TreePos->CurrentValue= TreePos->Child_One->CurrentValue
        / TreePos->Child_Two->CurrentValue;
      break;
    case MTK_POWER :
      TreePos->CurrentValue=pow(TreePos->Child_One->CurrentValue
                                ,TreePos->Child_Two->CurrentValue);
      break;
    case MTK_COMBINATION :
      TreePos->CurrentValue=combination(TreePos->Child_One->CurrentValue
                                        ,TreePos->Child_Two->CurrentValue);
      break;
    case MTK_PERMUTATION :
      TreePos->CurrentValue=permutation(TreePos->Child_One->CurrentValue
                                        ,TreePos->Child_Two->CurrentValue);
      break;
    case MTK_BESSJ :
      TreePos->CurrentValue=ansi_jn( (int)(TreePos->Child_One->CurrentValue) , TreePos->Child_Two->CurrentValue );
       break;
    default :
      throw ERR_FUNCTION;
    }
  }
  switch (errno) {
  case EDOM :
    errno=EZERO;
    throw ERR_DOM;
  case ERANGE :
    errno=EZERO;
    throw ERR_RANGE;
  };
  TreePos->Brackets=1;
}

// Return the name associated with a variable
const char *
Equation::VarText(int num)
{
  variable *temp;
  temp = Variables.at(num - 1);
  if (temp)
    return temp->Name;
  else 
    return "";
}

// Return the value associated with a variable
double
Equation::VarValue(int num)
{
  variable *temp;

  char memnan[8];
  memset(memnan, 0, 8);
  memnan[6] = (char) 0xF8;
  memnan[7] = (char) 0x7F;
  double *mynan = (double *)(memnan);
  
  temp = Variables.at(num - 1);
  if (temp)
    return temp->Value;
  else 
    return *mynan;
}

// Differentiate this equation and store in Differential
void 
Equation::Differentiate(Equation *Differential, int Var, Parser *parser)
{
  variable *PointerFrom, *PointerTo;
  DiffVar=Var;
  delete (Differential->Root);
  DiffFromNode(Root,&(Differential->Root));
  if (Differential->Root!=NULL) (Differential->Root)->Parent=NULL;
  PointerFrom = Variables.first();
  while ( PointerFrom != NULL ) {     
    PointerTo = new variable;
    PointerTo->Value = PointerFrom->Value;
    Differential->Variables.append(PointerTo);
    strcpy(PointerTo->Name,PointerFrom->Name);
    PointerFrom = Variables.next();
  };
  parser->TreeToText(Differential);
}

// Differentiate from this node
void 
Equation::DiffFromNode(Node *Start, Node **Hook)
{
  Node *Temp;

   switch (Start->Type)
   {
   case MTK_ADD : case MTK_SUBTRACT :
      *Hook=new Node;
      InitNode((*Hook),Start->Type);
      DiffFromNode(Start->Child_One,&((*Hook)->Child_One));
      ((*Hook)->Child_One)->Parent=(*Hook);
      DiffFromNode(Start->Child_Two,&((*Hook)->Child_Two));
      ((*Hook)->Child_Two)->Parent=(*Hook);
      break;
   case MTK_MULTIPLY :
      *Hook=new Node;
      InitNode((*Hook),MTK_ADD);
      AddBelowOne(*Hook,MTK_MULTIPLY);
      ((*Hook)->Child_One)->Child_One=new Node;
      Copy(Start->Child_One,((*Hook)->Child_One)->Child_One);
      (((*Hook)->Child_One)->Child_One)->Parent=(*Hook)->Child_One;
      DiffFromNode(Start->Child_Two,&(((*Hook)->Child_One)->Child_Two));
      (((*Hook)->Child_One)->Child_Two)->Parent=(*Hook)->Child_One;
      AddBelowTwo(*Hook,MTK_MULTIPLY);
      ((*Hook)->Child_Two)->Child_One=new Node;
      Copy(Start->Child_Two,((*Hook)->Child_Two)->Child_One);
      (((*Hook)->Child_Two)->Child_One)->Parent=(*Hook)->Child_Two;
      DiffFromNode(Start->Child_One,&(((*Hook)->Child_Two)->Child_Two));
      (((*Hook)->Child_Two)->Child_Two)->Parent=(*Hook)->Child_Two;
      break;
   case MTK_DIVIDE :
      *Hook=new Node;
      InitNode((*Hook),MTK_SUBTRACT);
      AddBelowOne(*Hook,MTK_DIVIDE);
      ((*Hook)->Child_One)->Child_Two=new Node;
      Copy(Start->Child_Two,((*Hook)->Child_One)->Child_Two);
      (((*Hook)->Child_One)->Child_Two)->Parent=(*Hook)->Child_One;
      DiffFromNode(Start->Child_One,&(((*Hook)->Child_One)->Child_One));
      (((*Hook)->Child_One)->Child_One)->Parent=(*Hook)->Child_One;
      AddBelowTwo(*Hook,MTK_DIVIDE);
      AddBelowOne((*Hook)->Child_Two,MTK_MULTIPLY);
      AddBelowTwo((*Hook)->Child_Two,MTK_POWER);
      (((*Hook)->Child_Two)->Child_One)->Child_One=new Node;
      Copy(Start->Child_One,(((*Hook)->Child_Two)->Child_One)->Child_One);
      ((((*Hook)->Child_Two)->Child_One)->Child_One)->Parent=((*Hook)->Child_Two)->Child_One;
      DiffFromNode(Start->Child_Two,&((((*Hook)->Child_Two)->Child_One)->Child_Two));
      ((((*Hook)->Child_Two)->Child_One)->Child_Two)->Parent=((*Hook)->Child_Two)->Child_One;
      (((*Hook)->Child_Two)->Child_Two)->Child_One=new Node;
      Copy(Start->Child_Two,(((*Hook)->Child_Two)->Child_Two)->Child_One);
      ((((*Hook)->Child_Two)->Child_Two)->Child_One)->Parent=((*Hook)->Child_Two)->Child_Two;
      AddBelowTwo(((*Hook)->Child_Two)->Child_Two,MTK_NUMBER,2);
      break;
   case MTK_POWER :
     *Hook = new Node;
     InitNode(*Hook, MTK_ADD);
     AddBelowOne(*Hook,MTK_MULTIPLY);
     AddBelowOne((*Hook)->Child_One,MTK_MULTIPLY);
     AddBelowTwo(*Hook,MTK_MULTIPLY);
     AddBelowOne((*Hook)->Child_Two,MTK_MULTIPLY);
     
     (*Hook)->Child_One->Child_One->Child_One = new Node;
     (*Hook)->Child_One->Child_One->Child_One->Parent = 
       (*Hook)->Child_One->Child_One;
     Copy(Start, (*Hook)->Child_One->Child_One->Child_One );
     AddBelowTwo((*Hook)->Child_One->Child_One, MTK_LN);
     (*Hook)->Child_One->Child_One->Child_Two->Child_One = new Node;
     (*Hook)->Child_One->Child_One->Child_Two->Child_One->Parent = 
       (*Hook)->Child_One->Child_One->Child_Two;
     Copy(Start->Child_One, (*Hook)->Child_One->Child_One->Child_Two->Child_One);
     DiffFromNode(Start->Child_Two,&((*Hook)->Child_One->Child_Two));     
     (*Hook)->Child_One->Child_Two->Parent
                = (*Hook)->Child_One;

     (*Hook)->Child_Two->Child_One->Child_One = new Node;
     (*Hook)->Child_Two->Child_One->Child_One->Parent = 
       (*Hook)->Child_Two->Child_One;
     Copy(Start, (*Hook)->Child_Two->Child_One->Child_One );
     Temp = new Node;
     InitNode(Temp, MTK_SUBTRACT);
     Temp->Parent = (*Hook)->Child_Two->Child_One->Child_One;
     Temp->Child_One = (*Hook)->Child_Two->Child_One->Child_One->Child_Two;
     Temp->Parent->Child_Two = Temp;
     Temp->Child_One->Parent = Temp;
     AddBelowTwo(Temp, MTK_NUMBER, 1);
     (*Hook)->Child_Two->Child_One->Child_Two = new Node;
     (*Hook)->Child_Two->Child_One->Child_Two->Parent = 
       (*Hook)->Child_Two->Child_One;
     Copy(Start->Child_Two, (*Hook)->Child_Two->Child_One->Child_Two);
     DiffFromNode(Start->Child_One,&((*Hook)->Child_Two->Child_Two));     
     (*Hook)->Child_Two->Child_Two->Parent = (*Hook)->Child_Two;
     break;
   case MTK_SIN :
      *Hook=new Node;
      InitNode((*Hook),MTK_MULTIPLY);
      AddBelowOne(*Hook,MTK_COS);
      ((*Hook)->Child_One)->Child_One=new Node;
      Copy(Start->Child_One,((*Hook)->Child_One)->Child_One);
      (((*Hook)->Child_One)->Child_One)->Parent=(*Hook)->Child_One;
      DiffFromNode(Start->Child_One,&((*Hook)->Child_Two));
      ((*Hook)->Child_Two)->Parent=(*Hook);
      break;
   case MTK_COS :
      (*Hook)=new Node;
      InitNode((*Hook),MTK_MULTIPLY);
      AddBelowOne((*Hook),MTK_MULTIPLY);
      AddBelowOne((*Hook)->Child_One,MTK_NUMBER,-1);
      AddBelowTwo((*Hook)->Child_One,MTK_SIN);
      (((*Hook)->Child_One)->Child_Two)->Child_One=new Node;
      Copy(Start->Child_One,(((*Hook)->Child_One)->Child_Two)->Child_One);
      ((((*Hook)->Child_One)->Child_Two)->Child_One)->Parent=((*Hook)->Child_One)->Child_Two;
      DiffFromNode(Start->Child_One,&((*Hook)->Child_Two));
      ((*Hook)->Child_Two)->Parent=(*Hook);
      break;
   case MTK_TAN :
      (*Hook)=new Node;
      InitNode((*Hook),MTK_MULTIPLY);
      AddBelowOne((*Hook),MTK_POWER);
      AddBelowOne((*Hook)->Child_One,MTK_SEC);
      AddBelowTwo((*Hook)->Child_One,MTK_NUMBER,2);
      (((*Hook)->Child_One)->Child_One)->Child_One=new Node;
      Copy(Start->Child_One,(((*Hook)->Child_One)->Child_One)->Child_One);
      ((((*Hook)->Child_One)->Child_One)->Child_One)->Parent=((*Hook)->Child_One)->Child_One;
      DiffFromNode(Start->Child_One,&((*Hook)->Child_Two));
      ((*Hook)->Child_Two)->Parent=(*Hook);
      break;
   case MTK_SEC :
      *Hook=new Node;
      InitNode((*Hook),MTK_MULTIPLY);
      AddBelowOne(*Hook,MTK_MULTIPLY);
      Temp = (*Hook)->Child_One;
      AddBelowOne(Temp,MTK_SEC);
      Temp->Child_One->Child_One=new Node;
      Copy(Start->Child_One,Temp->Child_One->Child_One);
      Temp->Child_One->Child_One->Parent=Temp->Child_One;
      AddBelowTwo(Temp,MTK_TAN);
      Temp->Child_Two->Child_One=new Node;
      Copy(Start->Child_One,Temp->Child_Two->Child_One);
      Temp->Child_Two->Child_One->Parent=Temp->Child_Two;
      DiffFromNode(Start->Child_One,&((*Hook)->Child_Two));
      ((*Hook)->Child_Two)->Parent=(*Hook);
      break;
   case MTK_COSEC :
      *Hook=new Node;
      InitNode((*Hook),MTK_MULTIPLY);
      AddBelowOne(*Hook,MTK_MULTIPLY);
      Temp = (*Hook)->Child_One;
      AddBelowOne(Temp,MTK_MULTIPLY);
      Temp = Temp->Child_One;
      AddBelowOne(Temp,MTK_NUMBER,-1);
      AddBelowTwo(Temp,MTK_COSEC);
      Temp->Child_Two->Child_One=new Node;
      Copy(Start->Child_One,Temp->Child_Two->Child_One);
      Temp->Child_Two->Child_One->Parent=Temp->Child_Two;
      Temp = Temp->Parent;
      AddBelowTwo(Temp,MTK_COT);
      Temp->Child_Two->Child_One=new Node;
      Copy(Start->Child_One,Temp->Child_Two->Child_One);
      Temp->Child_Two->Child_One->Parent=Temp->Child_Two;
      DiffFromNode(Start->Child_One,&((*Hook)->Child_Two));
      ((*Hook)->Child_Two)->Parent=(*Hook);
      break;
   case MTK_COT :
      (*Hook)=new Node;
      InitNode((*Hook),MTK_MULTIPLY);
      AddBelowOne(*Hook,MTK_MULTIPLY);
      Temp = (*Hook)->Child_One;
      AddBelowOne(Temp,MTK_NUMBER,-1);
      AddBelowTwo(Temp,MTK_POWER);
      Temp = Temp->Child_Two;
      AddBelowOne(Temp,MTK_COSEC);
      AddBelowTwo(Temp,MTK_NUMBER,2);
      Temp->Child_One->Child_One=new Node;
      Copy(Start->Child_One,Temp->Child_One->Child_One);
      Temp->Child_One->Child_One->Parent=Temp->Child_One;
      DiffFromNode(Start->Child_One,&((*Hook)->Child_Two));
      ((*Hook)->Child_Two)->Parent=(*Hook);
      break;
   case MTK_ARCSIN :
      (*Hook)=new Node;
      InitNode((*Hook),MTK_DIVIDE);
      DiffFromNode(Start->Child_One,&((*Hook)->Child_One));
      ((*Hook)->Child_One)->Parent=(*Hook);
      AddBelowTwo((*Hook),MTK_POWER);
      AddBelowOne((*Hook)->Child_Two,MTK_SUBTRACT);
      AddBelowOne(((*Hook)->Child_Two)->Child_One,MTK_NUMBER,1);
      AddBelowTwo(((*Hook)->Child_Two)->Child_One,MTK_POWER);
      ((((*Hook)->Child_Two)->Child_One)->Child_Two)->Child_One=new Node;
      Copy(Start->Child_One,((((*Hook)->Child_Two)->Child_One)->Child_Two)->Child_One);
      (((((*Hook)->Child_Two)->Child_One)->Child_Two)->Child_One)->Parent
                =(((*Hook)->Child_Two)->Child_One)->Child_Two;
      AddBelowTwo((((*Hook)->Child_Two)->Child_One)->Child_Two,MTK_NUMBER,2);
      AddBelowTwo((*Hook)->Child_Two,MTK_NUMBER,0.5);
      break;
   case MTK_ARCCOS :
      (*Hook)=new Node;
      InitNode((*Hook),MTK_DIVIDE);
      AddBelowOne((*Hook),MTK_MULTIPLY);
      AddBelowOne((*Hook)->Child_One,MTK_NUMBER,-1);
      DiffFromNode(Start->Child_One,&(((*Hook)->Child_One)->Child_Two));
      (((*Hook)->Child_One)->Child_Two)->Parent=(*Hook)->Child_One;
      AddBelowTwo((*Hook),MTK_POWER);
      AddBelowOne((*Hook)->Child_Two,MTK_SUBTRACT);
      AddBelowOne(((*Hook)->Child_Two)->Child_One,MTK_NUMBER,1);
      AddBelowTwo(((*Hook)->Child_Two)->Child_One,MTK_POWER);
      ((((*Hook)->Child_Two)->Child_One)->Child_Two)->Child_One=new Node;
      Copy(Start->Child_One,((((*Hook)->Child_Two)->Child_One)->Child_Two)->Child_One);
      (((((*Hook)->Child_Two)->Child_One)->Child_Two)->Child_One)->Parent
                =(((*Hook)->Child_Two)->Child_One)->Child_Two;
      AddBelowTwo((((*Hook)->Child_Two)->Child_One)->Child_Two,MTK_NUMBER,2);
      AddBelowTwo((*Hook)->Child_Two,MTK_NUMBER,0.5);
      break;
   case MTK_ARCTAN :
      (*Hook)=new Node;
      InitNode((*Hook),MTK_DIVIDE);
      DiffFromNode(Start->Child_One,&((*Hook)->Child_One));
      ((*Hook)->Child_One)->Parent=(*Hook);
      AddBelowTwo((*Hook),MTK_ADD);
      AddBelowOne((*Hook)->Child_Two,MTK_NUMBER,1);
      AddBelowTwo((*Hook)->Child_Two,MTK_POWER);
      (((*Hook)->Child_Two)->Child_Two)->Child_One=new Node;
      Copy(Start->Child_One,(((*Hook)->Child_Two)->Child_Two)->Child_One);
      ((((*Hook)->Child_Two)->Child_Two)->Child_One)->Parent
                =((*Hook)->Child_Two)->Child_Two;
      AddBelowTwo(((*Hook)->Child_Two)->Child_Two,MTK_NUMBER,2);
      break;
   case MTK_SINH :
      (*Hook)=new Node;
      InitNode((*Hook),MTK_MULTIPLY);
      AddBelowOne((*Hook),MTK_COSH);
      ((*Hook)->Child_One)->Child_One=new Node;
      Copy(Start->Child_One,((*Hook)->Child_One)->Child_One);
      (((*Hook)->Child_One)->Child_One)->Parent=(*Hook)->Child_One;
      DiffFromNode(Start->Child_One,&((*Hook)->Child_Two));
      ((*Hook)->Child_Two)->Parent=(*Hook);
      break;
   case MTK_COSH :
      (*Hook)=new Node;
      InitNode((*Hook),MTK_MULTIPLY);
      AddBelowOne((*Hook),MTK_SINH);
      ((*Hook)->Child_One)->Child_One=new Node;
      Copy(Start->Child_One,((*Hook)->Child_One)->Child_One);
      (((*Hook)->Child_One)->Child_One)->Parent=(*Hook)->Child_One;
      DiffFromNode(Start->Child_One,&((*Hook)->Child_Two));
      ((*Hook)->Child_Two)->Parent=(*Hook);
      break;
   case MTK_TANH :
      (*Hook)=new Node;
      InitNode((*Hook),MTK_DIVIDE);
      DiffFromNode(Start->Child_One,&((*Hook)->Child_One));
      ((*Hook)->Child_One)->Parent=(*Hook);
      AddBelowTwo((*Hook),MTK_POWER);
      AddBelowOne((*Hook)->Child_Two,MTK_COSH);
      AddBelowTwo((*Hook)->Child_Two,MTK_NUMBER,2);
      (((*Hook)->Child_Two)->Child_One)->Child_One=new Node;
      Copy(Start->Child_One,(((*Hook)->Child_Two)->Child_One)->Child_One);
      ((((*Hook)->Child_Two)->Child_One)->Child_One)->Parent=((*Hook)->Child_Two)->Child_One;
      break;
   case MTK_ARSINH :
      (*Hook)=new Node;
      InitNode(*Hook, MTK_DIVIDE);
      DiffFromNode(Start->Child_One,&((*Hook)->Child_One));
      ((*Hook)->Child_One)->Parent=(*Hook);
      AddBelowTwo((*Hook),MTK_POWER);
      AddBelowOne((*Hook)->Child_Two,MTK_ADD);
      AddBelowOne(((*Hook)->Child_Two)->Child_One,MTK_NUMBER,1);
      AddBelowTwo(((*Hook)->Child_Two)->Child_One,MTK_POWER);
      ((((*Hook)->Child_Two)->Child_One)->Child_Two)->Child_One=new Node;
      Copy(Start->Child_One,((((*Hook)->Child_Two)->Child_One)->Child_Two)->Child_One);
      (((((*Hook)->Child_Two)->Child_One)->Child_Two)->Child_One)->Parent
                =(((*Hook)->Child_Two)->Child_One)->Child_Two;
      AddBelowTwo((((*Hook)->Child_Two)->Child_One)->Child_Two,MTK_NUMBER,2);
      AddBelowTwo((*Hook)->Child_Two,MTK_NUMBER,0.5);
      break;
   case MTK_ARCOSH :
      (*Hook)=new Node;
      InitNode(*Hook, MTK_DIVIDE);
      DiffFromNode(Start->Child_One,&((*Hook)->Child_One));
      ((*Hook)->Child_One)->Parent=(*Hook);
      AddBelowTwo((*Hook),MTK_POWER);
      AddBelowOne((*Hook)->Child_Two,MTK_SUBTRACT);
      AddBelowTwo(((*Hook)->Child_Two)->Child_One,MTK_NUMBER,1);
      AddBelowOne(((*Hook)->Child_Two)->Child_One,MTK_POWER);
      ((((*Hook)->Child_Two)->Child_One)->Child_One)->Child_One=new Node;
      Copy(Start->Child_One,((((*Hook)->Child_Two)->Child_One)->Child_One)->Child_One);
      (((((*Hook)->Child_Two)->Child_One)->Child_One)->Child_One)->Parent
                =(((*Hook)->Child_Two)->Child_One)->Child_One;
      AddBelowTwo((((*Hook)->Child_Two)->Child_One)->Child_One,MTK_NUMBER,2);
      AddBelowTwo((*Hook)->Child_Two,MTK_NUMBER,0.5);
      break;
   case MTK_ARTANH :
      (*Hook)=new Node;
      InitNode(*Hook, MTK_DIVIDE);
      DiffFromNode(Start->Child_One,&((*Hook)->Child_One));
      ((*Hook)->Child_One)->Parent=(*Hook);
      AddBelowTwo((*Hook),MTK_SUBTRACT);
      AddBelowOne((*Hook)->Child_Two,MTK_NUMBER,1);
      AddBelowTwo((*Hook)->Child_Two,MTK_POWER);
      (((*Hook)->Child_Two)->Child_Two)->Child_One=new Node;
      Copy(Start->Child_One,(((*Hook)->Child_Two)->Child_Two)->Child_One);
      ((((*Hook)->Child_Two)->Child_Two)->Child_One)->Parent
                =((*Hook)->Child_Two)->Child_Two;
      AddBelowTwo(((*Hook)->Child_Two)->Child_Two,MTK_NUMBER,2);
      break;
   case MTK_LN :
      (*Hook)=new Node;
      InitNode(*Hook, MTK_DIVIDE);
      (*Hook)->Child_Two=new Node;
      Copy(Start->Child_One,(*Hook)->Child_Two);
      ((*Hook)->Child_Two)->Parent=(*Hook);
      DiffFromNode(Start->Child_One,&((*Hook)->Child_One));
      ((*Hook)->Child_One)->Parent=(*Hook);
      break;
   case MTK_LOG :
      (*Hook)=new Node;
      InitNode(*Hook, MTK_DIVIDE);
      AddBelowTwo(*Hook,MTK_MULTIPLY);
      AddBelowTwo((*Hook)->Child_Two,MTK_LN);
      AddBelowOne(((*Hook)->Child_Two)->Child_Two,MTK_NUMBER,10);
      ((*Hook)->Child_Two)->Child_One=new Node;
      Copy(Start->Child_One,((*Hook)->Child_Two)->Child_One);
      (((*Hook)->Child_Two)->Child_One)->Parent=(*Hook)->Child_Two;
      DiffFromNode(Start->Child_One,&((*Hook)->Child_One));
      ((*Hook)->Child_One)->Parent=*Hook;
      break;
   case MTK_EXP :
      (*Hook)=new Node;
      InitNode(*Hook, MTK_MULTIPLY);
      (*Hook)->Child_One=new Node;
      Copy(Start,(*Hook)->Child_One);
      ((*Hook)->Child_One)->Parent=(*Hook);
      DiffFromNode(Start->Child_One,&((*Hook)->Child_Two));
      ((*Hook)->Child_Two)->Parent=(*Hook);
      break;
   case MTK_NUMBER : case MTK_PI_SYM :
      (*Hook)=new Node;
      InitNode(*Hook, MTK_NUMBER);
      (*Hook)->Value.Number = 0;
      (*Hook)->Child_One=NULL;
      (*Hook)->Child_Two=NULL;
      break;
   case MTK_VARIABLE :
      (*Hook)=new Node;
      InitNode(*Hook, MTK_NUMBER);
      (*Hook)->Child_One=NULL;
      (*Hook)->Child_Two=NULL;
      if ((Start->Value).WhichVar==DiffVar)
      {
         ((*Hook)->Value).Number=1;
      }
      else
      {
         ((*Hook)->Value).Number=0;
      };
      break;
   default :
      throw ERR_DIFF;
      /*(*Hook)=new Node;
      InitNode(*Hook, MTK_NUMBER);
      (*Hook)->Child_One=NULL;
      (*Hook)->Child_Two=NULL;
      break;*/
   };

}

// Numerically integrate equation
double
Equation::NumInt( range Range, int var_num )
{
  variable *var;
  double n,Ans;
  int i;

  var = Variables.at( var_num - 1 );

  if (!var) throw ERR_VARIABLES;

  n = floor( ( Range.To - Range.From ) / ( Range.Step * 2 ) );
  if ( n==0 ) n=1;
  Range.Step = ( Range.To - Range.From ) / ( n * 2 );
  var->Value = Range.From;
  Ans = Evaluate();
  for ( i=1 ; i<n ; i++ ) {
    var->Value = Range.From + ( 2*i-1 ) * Range.Step;
    Ans = Ans + 4 * Evaluate();
    var->Value = var->Value + Range.Step;
    Ans = Ans + 2 * Evaluate();
  };
  var->Value = var->Value + Range.Step;
  Ans = Ans + 4 * Evaluate();
  var->Value = Range.To;
  Ans = Ans + Evaluate();
  Ans = Ans * Range.Step / 3;
  return Ans;
}

// Find a root of equation in interval by bisection
double 
Equation::FindRoot( range Range, int var_num )
{
  double Mid,Below,Above,AboveVal,BelowVal,MidVal,dx;
  variable *var;
  int i;

  Below = Range.From;
  Above = Range.To;
  var = Variables.at( var_num - 1 );
  var->Value = Below;
  BelowVal = Evaluate();
  if ( fabs(BelowVal) < Range.Accuracy ) return Below;
  var->Value = Above;
  AboveVal = Evaluate();
  if ( fabs(AboveVal) < Range.Accuracy ) return Above;

  // If no root sign change range then try subdividing range
  if ( AboveVal*BelowVal > 0 ) {
    dx = (Range.To - Range.From)/50;
    Below = Range.From;
    Above = Below + dx;
    for (i=0;i<50;i++) {
      var->Value = Below;
      BelowVal = Evaluate();
      var->Value = Above;
      AboveVal = Evaluate();
      if ( AboveVal*BelowVal < 0 ) break;
      Below += dx;
      Above += dx;
    }

    // No success so try expanding
    if ( AboveVal*BelowVal > 0 ) {
      Below = Range.From;
      Above = Range.To;
      var->Value = Below;
      BelowVal = Evaluate();
      var->Value = Above;
      AboveVal = Evaluate();
      for (i=0; i<50;i++) {
        if ( AboveVal*BelowVal < 0 ) break;      
        if ( fabs(BelowVal) < fabs(AboveVal) ) {
          Below += 1.6*(Below-Above);
          var->Value = Below;
          BelowVal = Evaluate();
        }
        else {
          Above += 1.6*(Above-Below);
          var->Value = Above;
          AboveVal = Evaluate();
        }
      }
      // Still no luck - throw an error
      if ( AboveVal*BelowVal > 0 ) {
        throw ERR_NO_ROOT;
      }
    }
  }

  do {
    Mid = ( Above + Below ) / 2;
    var->Value = Mid;
    MidVal = Evaluate();
    if ( MidVal*BelowVal < 0 ) {
      Above = Mid;
      AboveVal = MidVal;
    }
    else {
      Below = Mid;
      BelowVal = MidVal;
    };
  }
  while ( fabs(MidVal) > Range.Accuracy );
  return Mid;
}

// Solve ODE using fourth order Runge-Kutta
double 
Equation::RK4( range Range, int var_num1, int var_num2 )
{
  variable *var1=NULL;
  variable *var2=NULL;
  double n,x,ans,hh;
  double k1,k2,k3,k4;
  int i;
  if (var_num1) var1 = Variables.at( var_num1 - 1 );
  if (var_num2) var2 = Variables.at( var_num2 - 1 );

  hh = Range.Step/2;
  n = floor( (Range.To-Range.From) / Range.Step );
  if ( n==0 ) n=1;
  n = fabs( n );
  Range.Step = (Range.To-Range.From) / n;
  x = Range.From;
  ans = Range.Accuracy;
  for ( i=1; i <= (int) n; i++ ) {                
    if ( var1 ) var1->Value = x;
    if ( var2 ) var2->Value = ans;
    k1 = Range.Step * Evaluate();
    if ( var1 ) var1->Value = x + hh;
    if ( var2 ) var2->Value = ans + k1 / 2;
    k2 = Range.Step * Evaluate();
    if ( var2 ) var2->Value = ans + k2 / 2;
    k3 = Range.Step * Evaluate();
    if ( var1 ) var1->Value = x + Range.Step;
    if ( var2 ) var2->Value = ans + k3;
    k4 = Range.Step * Evaluate();
    ans = ans + ( k1 / 2 + k2 + k3 + k4 / 2 ) / 3;
    x = x + Range.Step;
  };
  return ans;
}

// Calculate first n terms of power series
void 
Equation::PowerSeries(Equation *Series, int var, int n, Parser *parser )
{
  Equation *Differential=NULL;
  Equation *Temp=NULL;
  variable *PointerFrom;
  variable *PointerTo;
  variable *AboutVar=NULL;
  int i;
  double AboutValue;

  AboutVar = Variables.at( var - 1 );
  if (!AboutVar) {
    throw ERR_VARIABLES;
  }

  AboutValue = AboutVar->Value;
  
  // Copy variables to series
  PointerFrom = Variables.first();
  while ( PointerFrom != NULL ) {     
    PointerTo = new variable;
    PointerTo->Value = PointerFrom->Value;
    Series->Variables.append(PointerTo);
    strcpy(PointerTo->Name,PointerFrom->Name);
    PointerFrom = Variables.next();
  };

  Series->Radians = Radians;
  if (Series->Root) {
    delete Series->Root;
    Series->Root = NULL;
  }
  Series->AddPowerSeriesTerm(var, Root, 0, AboutValue,parser);
  Temp = this;

  for (i=1;i<n;i++) {
    DiffVar=var;
    Differential = new Equation();
    delete (Differential->Root);
    Differential->Root = NULL;
    DiffFromNode(Temp->Root,&(Differential->Root));
    if (Differential->Root!=NULL) Differential->Root->Parent=NULL;
    PointerFrom = Variables.first();
    while ( PointerFrom != NULL ) {     
      PointerTo = new variable;
      PointerTo->Value = PointerFrom->Value;
      Differential->Variables.append(PointerTo);
      strcpy(PointerTo->Name,PointerFrom->Name);
      PointerFrom = Variables.next();
    };
    parser->TidyUpEqn(Differential);
    Series->AddPowerSeriesTerm(var, Differential->Root, i, AboutValue, parser);
    if ( Temp != this ) delete Temp;
    Temp = Differential;
  }
  if ( Temp != this ) delete Temp;
  parser->TreeToText(Series);
}

// Add power series term to equation
void
Equation::AddPowerSeriesTerm(int var, Node *coeff, double power, double about, Parser *parser)
{
  Node *temp;
  Node *temp_root;
  NumVar tmpnum;

  if (Root != NULL) {
    temp = new Node;
    temp->Parent = NULL;
    temp->Child_One = Root;
    Root->Parent = temp;
    InitNode(temp,MTK_ADD);
    Root = temp;

    temp->Child_Two = new Node;
    temp->Child_Two->Parent = temp;
    temp = temp->Child_Two;
  }
  else {
    temp = new Node;
    Root = temp;
    temp->Parent = NULL;
  }
  
  InitNode(temp,MTK_MULTIPLY);
  temp->Child_One = new Node;
  temp->Child_One->Parent = temp;
  temp = temp->Child_One;
  InitNode(temp,MTK_DIVIDE);

  temp_root = new Node;
  temp->Child_One = temp_root;
  temp_root->Parent = temp;
  Copy(coeff,temp_root);

  temp->Child_Two = new Node;
  temp->Child_Two->Parent = temp;
  tmpnum.Number = factorial(power);
  InitNode(temp->Child_Two,MTK_NUMBER,tmpnum);
  temp->Child_Two->Child_One = NULL;
  temp->Child_Two->Child_Two = NULL;
  temp = temp->Parent;
    
  temp->Child_Two = new Node;
  temp->Child_Two->Parent = temp;
  temp = temp->Child_Two;
  InitNode(temp,MTK_POWER);

  temp->Child_One = new Node;
  temp->Child_One->Parent = temp;
  temp = temp->Child_One;
  InitNode(temp,MTK_SUBTRACT);

  temp->Child_One = new Node;
  temp->Child_One->Parent = temp;
  temp = temp->Child_One;
  tmpnum.WhichVar = var;
  InitNode(temp,MTK_VARIABLE,tmpnum);
  temp->Child_One = NULL;
  temp->Child_Two = NULL;
  temp = temp->Parent;

  temp->Child_Two = new Node;
  temp->Child_Two->Parent = temp;
  temp = temp->Child_Two;
  tmpnum.Number = about;
  InitNode(temp,MTK_NUMBER,tmpnum);
  temp->Child_One = NULL;
  temp->Child_Two = NULL;

  temp = temp->Parent->Parent;
  temp->Child_Two = new Node;
  temp->Child_Two->Parent = temp;
  temp = temp->Child_Two;
  tmpnum.Number = power;
  InitNode(temp,MTK_NUMBER,tmpnum);
  temp->Child_One = NULL;
  temp->Child_Two = NULL;

  temp = temp_root;
  Reset(0,temp_root);
  parser->TidyUpEqn(this);
  Reset(0);
  while (temp_root->Brackets==0) {
    while (temp->Child_One!=NULL) {
      if (temp->Child_One->Brackets==1) break;
      temp=temp->Child_One;
    };
    if (temp->Type == MTK_VARIABLE) {
      if (temp->Value.WhichVar == var) {
        temp->Type = MTK_NUMBER;
        temp->evalfn = func_details[MTK_NUMBER].func_ptr;
        temp->Value.Number = Variables.at(var - 1)->Value;
      }
    }
    temp->Brackets=1;
    temp=temp->Parent;
    if (temp==NULL) break;
    if (temp->Child_Two != NULL) {
      if ((temp->Child_Two)->Brackets==0) {
        temp=temp->Child_Two;
      }
      else {
        temp->Brackets = 1;
      }
    }
    else {
      temp->Brackets = 1;
    };
  };
  Reset(2);
  Reset(0);
}
