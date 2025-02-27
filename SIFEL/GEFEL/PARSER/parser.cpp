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

// Parser class

#define PARSER_CPP 1

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include "moremath.h"
#include "parser.h"
#include "equ.h"
#include "func.h"
#include "functext.h"
#include "errors.h"

//Constructor for the Parser class
Parser::Parser()
{
  EquationTokenStart = NULL;
  EquationTokenEnd = NULL;
  EquationText = NULL;
  equation = NULL;
}

//Destructor for the Parser class
Parser::~Parser()
{
  while (EquationTokenStart!=EquationTokenEnd) {
    EquationTokenStart=EquationTokenStart->Next;
    delete EquationTokenStart->Last;
  };
  if (EquationTokenStart!=NULL) delete EquationTokenStart;
}

// Convert text to a tree and return error code
Equation *
Parser::TextToTree( char *text )
{
  try {
    equation = new Equation();

    EquationText = new char[strlen(text) + 1];
    strcpy(EquationText,text);
    if (strlen(EquationText) == 0) {
      throw ERR_NO_TEXT;
    };
    PassOne();
    ErrorCheck();
    PassTwo();
    RemoveBrackets();
    TidyUp();
    TreeToText();
  }
  catch (...) {
    if (equation) 
    {
      delete equation;
      equation = NULL;
    }
    if (EquationText) {
      delete [] EquationText;
      EquationText = NULL;
    }
    while (EquationTokenStart!=EquationTokenEnd) {
      EquationTokenStart=EquationTokenStart->Next;
      delete EquationTokenStart->Last;
    };
    if (EquationTokenStart!=NULL) {
      delete EquationTokenStart;
      EquationTokenStart = NULL;
      EquationTokenEnd = NULL;
    }
    throw;
  }

  if (equation->EquationText)
    delete [] equation->EquationText;
  equation->EquationText = EquationText;
  EquationText = NULL;
  while (EquationTokenStart!=EquationTokenEnd) {
    EquationTokenStart=EquationTokenStart->Next;
    delete EquationTokenStart->Last;
  };
  if (EquationTokenStart!=NULL) {
    delete EquationTokenStart;
    EquationTokenStart = NULL;
    EquationTokenEnd = NULL;
  }
  
  return equation;
}

// Convert tree to text
void
Parser::TreeToText( Equation *eqn )
{
  equation = eqn;

  try {
    TidyUp();
    TreeToText();
  }
  catch (...) {
    if (EquationText) {
      delete [] EquationText;
      EquationText = NULL;
    }
  }
  
  if (equation->EquationText != NULL) delete [] equation->EquationText;
  equation->EquationText = EquationText;
  EquationText = NULL;
}

void 
Parser::TidyUpEqn(Equation *eqn)
{
  equation = eqn;

  TidyUp();
}

// First pass of the parser
void 
Parser::PassOne()
{
  char * Pointer;
  Pointer = EquationText;
  char * ThisFunc = new char[10];
  ThisFunc[0] = 0;
  double * ThisDigit = new double;
  *ThisDigit = 0;
  bool * DoDigit = new bool;
  *DoDigit = false;
  int * DecPt = new int;
  *DecPt = 0;
  double * ThisExp = new double;
  *ThisExp = 0;
  bool * DoDecExp = new bool;
  *DoDecExp = false;
  bool ExpSgn = false;
  bool ReadExpSgn = false;
  int dummy = 0;
  
  equation->brackets = 0;
  
  while ( *Pointer != 0 ) {
    if ( isdigit(*Pointer) ) {
      if ( strlen(ThisFunc) != 0 ) EndFunc(ThisFunc);
      if (*DoDecExp == false)
        Digit( Pointer, ThisDigit, DoDigit, DecPt );
      else
      {
        ReadExpSgn = false;
        Digit( Pointer, ThisExp, DoDecExp, &dummy );
      }
    }
    if (tolower(*Pointer) == 'e')
    {
      if ((*DoDigit == true) && (*DoDecExp == false))
      {
        *DoDecExp = true;
        ExpSgn = false;
        ReadExpSgn = true;
        Pointer++;
        continue;
      }
    }
    if ( isalpha(*Pointer) ) {
      if ( *DoDigit == true ) {      //End of number
        if (ExpSgn == true) {
          *ThisExp = -(*ThisExp);
          ExpSgn = false;
        }
        *ThisDigit *= pow(10.0, *ThisExp);
        *ThisExp = 0.0;  *DoDecExp = false;
        EndNum( ThisDigit, DoDigit, DecPt );
      }
      Alpha( Pointer, ThisFunc );
    }
    if ( ispunct(*Pointer) ) {
      if ( ( *DoDigit == true ) && ( *Pointer != '.' ) && ( *DoDecExp == false)) {      //Not Decimal Pt, end of number
        if (ExpSgn == true) {
          *ThisExp = -(*ThisExp);
          ExpSgn = false;
        }
        *ThisDigit *= pow(10.0, *ThisExp);
        *ThisExp = 0.0;  *DoDecExp = false;
        EndNum( ThisDigit, DoDigit, DecPt );
      }
      if (( *DoDigit == true ) && ( *DoDecExp == true) && (ReadExpSgn == false)) {  //not sign of exponent, end of number
        if (ExpSgn == true) {
          *ThisExp = -(*ThisExp);
          ExpSgn = false;
        }
        *ThisDigit *= pow(10.0, *ThisExp);
        *ThisExp = 0.0;  *DoDecExp = false;
        EndNum( ThisDigit, DoDigit, DecPt );
      }
      if ( strlen(ThisFunc) != 0 ) EndFunc(ThisFunc);
      if ( *Pointer == '.' ) {     //Decimal Pt in number
        *DecPt = 1;
        if ( *DoDigit == false ) {
          *DoDigit = true;
          *ThisDigit = 0;
        }
      }
      else {
        if ((*Pointer == '-') && (*DoDecExp == true) && (ReadExpSgn == true))
        {
          ExpSgn = true;
          ReadExpSgn = false;
          Pointer++;
          continue;
        }
        if ((*Pointer == '+') && (*DoDecExp == true) && (ReadExpSgn == true))
        {
          ReadExpSgn = false;
          Pointer++;
          continue;
        }
        Symbol(Pointer);
      }
    }
    if (isspace(*Pointer)) {
      if ( *DoDigit == true ) 
      {
        if (ExpSgn == true) {
          *ThisExp = -(*ThisExp);
          ExpSgn = false;
        }
        *ThisDigit *= pow(10.0, *ThisExp);
        *ThisExp = 0.0;  *DoDecExp = false;
        EndNum( ThisDigit, DoDigit, DecPt );
      }
      if ( ThisFunc[0] != 0 ) EndFunc( ThisFunc );
    }
    Pointer++;
  }
  if ( *DoDigit == true ) 
  {
    if (ExpSgn == true) {
      *ThisExp = -(*ThisExp);
      ExpSgn = false;
    }
    *ThisDigit *= pow(10.0, *ThisExp);
    *ThisExp = 0.0;  *DoDecExp = false;
    EndNum( ThisDigit, DoDigit, DecPt );
  }
  if ( ThisFunc[0] != 0 ) EndFunc( ThisFunc );
  delete ThisDigit;
  delete DoDigit;
  delete DecPt;
  delete DoDecExp;
  delete ThisExp; 
  delete [] ThisFunc;  
}

// Add next digit to number
void 
Parser::Digit(char * Pointer, double * ThisDigit, bool * DoDigit, int * DecPt)
{
   if ( *DecPt == 0 ) {
      *ThisDigit = (*ThisDigit) * 10 + (*Pointer) - 48;
   }
   else {
     (*ThisDigit) = (*ThisDigit) + ((*Pointer) - 48) / pow( 10.0, (double)(*DecPt) );
     (*DecPt)++;
   }
   if ( *DoDigit == false ) *DoDigit = true;
}

// Add alphabetic character to this function string
void 
Parser::Alpha(const char * Pointer, char * ThisFunc)
{
  if ( strlen(ThisFunc) < 10 ) {
    strncat( ThisFunc, Pointer, 1 );
  }
}

// Check to see what a symbol is
void 
Parser::Symbol(char * Pointer)
{
  func_type WhichSymbol = MTK_NOFUNC;
  switch (*Pointer) {
  case '!' : 
    WhichSymbol=MTK_FACTORIAL;   // !
    break;
  case '(' : WhichSymbol=MTK_LBRACKET;    // (
    equation->brackets++;
    break;
  case ')' : WhichSymbol=MTK_RBRACKET;    // )
    equation->brackets--;
    break;
  case '*' : WhichSymbol=MTK_MULTIPLY;    // *
    break;
  case '+' : WhichSymbol=MTK_ADD;         // +
    break;
  case ',': WhichSymbol = MTK_COMMA;      // ,
    break;
  case '-' : WhichSymbol=MTK_SUBTRACT;    // -
    break;
  case '/' : WhichSymbol=MTK_DIVIDE;      // /
    break;
  case '^' : WhichSymbol=MTK_POWER;       // ^
    break;
  default :         // Anything else not recognised
    // Tidy up then throw an error
    throw ERR_SYMBOL;
  };
  if (WhichSymbol!=0) {
    AddToken(WhichSymbol);
    (EquationTokenEnd->Value).Number=0;
  };
}

// Check to see if function string is a function
void 
Parser::EndFunc(char * ThisFunc)
{
  func_type WhichFunc=MTK_NOFUNC, i;

  if (strcmp(ThisFunc,"C")==0) {
    WhichFunc=MTK_COMBINATION;
  }
  else {
    if (strcmp(ThisFunc,"P")==0) {
      WhichFunc=MTK_PERMUTATION;
    }
    else {
      if (strcmp(ThisFunc,"pi")==0) WhichFunc=MTK_PI_SYM;
      else {
        Alpha(" ",ThisFunc);
        i = MTK_LN;
        while ( WhichFunc==MTK_NOFUNC && i<=MTK_FRAC ) {
          if ( strcmp(ThisFunc,func_details[i].func_text)==0 ) WhichFunc = i;
          i++;
        }
      }
    }
  }

  if (WhichFunc != MTK_NOFUNC) {
    AddToken(WhichFunc);
    EquationTokenEnd->Value.Number=0;
  }
  else {
    variable *VarPoint;
    VarPoint = equation->Variables.first();
    AddToken(MTK_VARIABLE);
    ThisFunc[ strlen(ThisFunc) - 1 ] = 0;
    while (VarPoint!=NULL) {
      if (strcmp(ThisFunc,VarPoint->Name)==0) {
        EquationTokenEnd->Value.WhichVar=equation->Variables.at() + 1;
        ThisFunc[0]=0;
        return;
      }
      VarPoint=equation->Variables.next();
    }
    VarPoint = new variable;
    equation->Variables.append(VarPoint);
    EquationTokenEnd->Value.WhichVar=equation->Variables.count();
    strncpy( VarPoint->Name, ThisFunc, 12 );
    VarPoint->Value = 0;
  }
  ThisFunc[0] = 0;
}

// End of number reached - add token
void 
Parser::EndNum(double * ThisDigit, bool * DoDigit, int * DecPt)
{
    AddToken(MTK_NUMBER);
    (EquationTokenEnd->Value).Number=*ThisDigit;
    *ThisDigit=0;
    *DoDigit=false;
    *DecPt=0;
}

// Add a token to the token representation
void 
Parser::AddToken(func_type WhichType)
{
    if (EquationTokenEnd!=NULL)
    {
       EquationTokenEnd->Next=new Token;
       (EquationTokenEnd->Next)->Last=EquationTokenEnd;
       EquationTokenEnd=EquationTokenEnd->Next;
       EquationTokenEnd->Next=NULL;
       EquationTokenEnd->Type=WhichType;
    }
    else
    {
       EquationTokenStart=new Token;
       EquationTokenEnd=EquationTokenStart;
       EquationTokenStart->Last=NULL;
       EquationTokenStart->Next=NULL;
       EquationTokenStart->Type=WhichType;
    };
}



// Check for various errors
void 
Parser::ErrorCheck()
{
  Token *Pointer,*Temp;

  if (equation->brackets != 0) {
    throw ERR_BRACKET;
  }

  Pointer=EquationTokenStart;
  if (EquationTokenStart==EquationTokenEnd) {
    if ( EquationTokenStart==NULL ) {
      throw ERR_NO_TEXT;
    };
    if (((Pointer->Type)<MTK_NUMBER)||((Pointer->Type)>MTK_PI_SYM)) {
      throw ERR_FUNCTION;
    };
    return;
  };
  switch (Pointer->Type) {
  case MTK_ADD :     
    if (((((Pointer->Next)->Type)>=MTK_MULTIPLY)&&(((Pointer->Next)->Type)<=MTK_POWER)) || (Pointer->Next->Type == MTK_BESSJ)) {
      throw ERR_OPERATOR;
    };
    (Pointer->Next)->Last=NULL;
    Temp=Pointer;
    Pointer=Pointer->Next;
    delete Temp;
    EquationTokenStart=Pointer;
    break;
  case MTK_SUBTRACT :
    if (((((Pointer->Next)->Type)>=MTK_MULTIPLY)&&(((Pointer->Next)->Type)<=MTK_POWER)) || (Pointer->Next->Type == MTK_BESSJ)) {
      throw ERR_OPERATOR;
    };
    if ((Pointer->Next)->Type==MTK_NUMBER) {
      (((Pointer->Next)->Value).Number)=-(((Pointer->Next)->Value).Number);
      (Pointer->Next)->Last=NULL;
      Temp=Pointer;
      Pointer=Pointer->Next;
      delete Temp;
      EquationTokenStart=Pointer;
    }
    else {
      Pointer->Type=MTK_NUMBER;
      (Pointer->Value).Number=-1;
      Temp=Pointer->Next;
      Pointer->Next = new Token;
      Pointer->Next->Last = Pointer;
      Pointer->Next->Next = Temp;
      Pointer->Next->Type = MTK_MULTIPLY;
      Temp->Last = Pointer->Next;
    };
    break;
  case MTK_MULTIPLY : case MTK_DIVIDE : case MTK_POWER :
    throw ERR_OPERATOR;
  case MTK_COMBINATION : case MTK_PERMUTATION : case MTK_FACTORIAL : case MTK_COMMA :
    throw ERR_FUNCTION;
  case MTK_RBRACKET :
    throw ERR_BRACKET;
  };
  while ((Pointer->Next!=EquationTokenEnd)&&(Pointer->Next!=NULL)) {
    Pointer=Pointer->Next;
    switch (Pointer->Type) {
    case MTK_MULTIPLY : case MTK_DIVIDE : case MTK_POWER : 
      if ( (Pointer->Last->Type != MTK_LBRACKET) && ((Pointer->Last->Type < MTK_FACTORIAL) || (Pointer->Last->Type > MTK_RBRACKET)) ) {
        throw ERR_OPERATOR;
      };
      break;
    case MTK_FACTORIAL : case MTK_COMBINATION : case MTK_PERMUTATION : case MTK_COMMA :
      if ( (Pointer->Last->Type < MTK_NUMBER) || (Pointer->Last->Type > MTK_RBRACKET) ) {
        if (Pointer->Last->Type != MTK_INT) {
          throw ERR_FUNCTION;
        }
      }
      break;
    case MTK_ADD : case MTK_SUBTRACT :
      if (((Pointer->Next->Type >= MTK_MULTIPLY) && (Pointer->Next->Type <= MTK_POWER)) || (Pointer->Next->Type == MTK_COMMA)) {
        throw ERR_OPERATOR;
      };
      if ( ( (Pointer->Last->Type < MTK_NUMBER) || (Pointer->Last->Type > MTK_RBRACKET) ) && (Pointer->Last->Type != MTK_FACTORIAL) ){
        if ( (Pointer->Next->Type == MTK_NUMBER) && (Pointer->Type == MTK_SUBTRACT) ) 
          Pointer->Next->Value.Number = - Pointer->Next->Value.Number;
        if ( (Pointer->Next->Type == MTK_NUMBER) || (Pointer->Type == MTK_ADD) ) {
          Pointer->Next->Last = Pointer->Last;
          Temp = Pointer;
          Pointer = Pointer->Next;
          delete Temp;
          Pointer->Last->Next = Pointer;
        }
        else {
          Pointer->Type=MTK_NUMBER;
          (Pointer->Value).Number=-1;
          Temp=Pointer->Next;
          Pointer->Next = new Token;
          Pointer->Next->Last = Pointer;
          Pointer->Next->Next = Temp;
          Pointer->Next->Type = MTK_MULTIPLY;
          Temp->Last = Pointer->Next;
        };
      };
      break;
    case MTK_RBRACKET :
      if ((((Pointer->Last)->Type)<MTK_LBRACKET)
          || (((Pointer->Last)->Type)>MTK_RBRACKET)) {
        throw ERR_BRACKET;
      };
      break;
    default :
      if ((Pointer->Last->Type > MTK_LBRACKET) && (Pointer->Last->Type < MTK_COMMA)) {
        Temp = Pointer->Last;
        Pointer->Last=new Token;
        Pointer->Last->Next=Pointer;
        Pointer=Pointer->Last;
        Pointer->Last=Temp;
        Temp->Next=Pointer;
        Pointer->Type=MTK_MULTIPLY;
        Pointer->Value.Number=0;
        Pointer=Pointer->Next;
      };
    };
  };
  Pointer=Pointer->Next;
  if (Pointer) {
    if ((((Pointer->Type)>=MTK_COMBINATION)&&((Pointer->Type)<=MTK_FRAC)) 
        || (Pointer->Type == MTK_COMMA) ) {
      throw ERR_FUNCTION;
    };
    if ((Pointer->Type)==MTK_LBRACKET) {
      throw ERR_BRACKET;
    };
    if (((Pointer->Type)<=MTK_POWER) && ((Pointer->Type)>=MTK_ADD)) {
      throw ERR_OPERATOR;
    };
    if (((Pointer->Type)>=MTK_NUMBER)&&((Pointer->Type)<=MTK_PI_SYM)) {
      if ((Pointer->Last->Type >= MTK_NUMBER) && (Pointer->Last->Type < MTK_COMMA)) {
        Temp=Pointer->Last;
        Pointer->Last=new Token;
        Pointer->Last->Next=Pointer;
        Pointer=Pointer->Last;
        Pointer->Last=Temp;
        Temp->Next=Pointer;
        Pointer->Type=MTK_MULTIPLY;
        Pointer->Value.Number=0;
        Pointer=Pointer->Next;
      }
    }
  }
}


// Second pass of equation parser
void 
Parser::PassTwo()
{
   Node *TreePos;
   Token *Pointer;

   Pointer = EquationTokenStart;
   TreePos = equation->Root;
   TreePos->Type=Pointer->Type;
   TreePos->evalfn = func_details[TreePos->Type].func_ptr;
   TreePos->Value = Pointer->Value;

   while (Pointer->Next!=NULL) {
     Pointer=Pointer->Next;

     if (Priority(TreePos->Type,Pointer->Type)==true) {
       if ((TreePos->Parent!=NULL)&&(Pointer->Type!=MTK_FACTORIAL)) {
         while (Priority((TreePos->Parent)->Type,Pointer->Type)==true) {
           TreePos=TreePos->Parent;
           if (TreePos->Parent==NULL) break;
         }
       }
       if (Pointer->Type==MTK_RBRACKET) {
         if (TreePos->Parent==NULL) {
           throw ERR_BRACKET;
         }
         TreePos=TreePos->Parent;
       }
       TreePos=equation->AddParent(TreePos);
     }
     else {
       if ((((TreePos->Type>=MTK_ADD)&&(TreePos->Type)<=MTK_POWER))||(Pointer->Type==MTK_LBRACKET) || (TreePos->Type == MTK_COMMA)) {
         if (Pointer->Type==MTK_LBRACKET) {
           if (TreePos->Child_One==NULL) {
             TreePos=equation->AddChildOne(TreePos);
           }
           else {
             TreePos=equation->AddChildTwo(TreePos);
           }
         }
         else {
           TreePos=equation->AddChildTwo(TreePos);
         }
       }
       else {
         while ((Priority(TreePos->Type,Pointer->Type)!=true)
                ||(TreePos->Child_One!=NULL)) {
           if (TreePos->Child_One==NULL) break;
           if (TreePos->Child_Two==NULL) {
             TreePos=TreePos->Child_One;
           }
           else {
             TreePos=TreePos->Child_Two;
           }
         }
         if (Priority(TreePos->Type,Pointer->Type)==true) {
           TreePos=equation->AddParent(TreePos);
         }
         else {
           TreePos=equation->AddChildOne(TreePos);
         }
       }
     }
     equation->InitNode(TreePos,Pointer->Type,Pointer->Value);
     TreePos->Child_Two=NULL;
   }
}

//Returns true if First is higher or equal priority than Second else returns false
bool 
Parser::Priority(func_type First, func_type Second)         
{
    if (Second==MTK_LBRACKET) return false;
    if (First==MTK_LBRACKET) return false;
    if (First==MTK_RBRACKET) return true;
    if (Second==MTK_RBRACKET) return true;
    if (Second == MTK_COMMA) return true;
    if (First == MTK_COMMA) return false;
    if ((First==MTK_ADD)||(First==MTK_SUBTRACT))
    {
       if ((Second==MTK_ADD)||(Second==MTK_SUBTRACT)) return true;
       return false;
    };
    if ((First==MTK_MULTIPLY)||(First==MTK_DIVIDE))
    {
       if ((Second>=MTK_ADD)&&(Second<=MTK_DIVIDE)) return true;
       return false;
    };
    if (Second==MTK_FACTORIAL) return true;
    if ((Second>=MTK_ADD)&&(Second<=MTK_DIVIDE)) return true;
    if (Second == MTK_BESSJ) return true;
    if ((First>=MTK_NUMBER)&&(First<=MTK_RBRACKET)) return true;
    return false;
}

// Remove unnecessary brackets from the equation
void
Parser::RemoveBrackets()
{
   Node *TreePos, *temp;

   TreePos = equation->Root;
   while (equation->Root->Brackets==0) {
     while (TreePos->Child_One!=NULL) {
       if ((TreePos->Child_One)->Brackets==1) break;
       TreePos=TreePos->Child_One;
     };
     if ((TreePos->Type==MTK_LBRACKET)||(TreePos->Type==MTK_RBRACKET)) {
       TreePos=equation->DeleteNode(TreePos);
       if (TreePos==NULL) break;
     }
     else {
       if (TreePos->Child_Two==NULL) {
         TreePos->Brackets=1;
         TreePos=TreePos->Parent;
       }
       else {
         if ((TreePos->Child_Two)->Brackets==1) {
           TreePos->Brackets=1;
           TreePos=TreePos->Parent;
         }
         else {
           TreePos=TreePos->Child_Two;
         }
       }
     }
   }
   if ((equation->Root->Type==MTK_LBRACKET)||(equation->Root->Type==MTK_RBRACKET)) {
     equation->Delete(equation->Root);
   }

   equation->Reset(0);
   TreePos = equation->Root;
   while (equation->Root->Brackets==0) {
     while (TreePos->Child_One!=NULL) {
       if ((TreePos->Child_One)->Brackets==1) break;
       TreePos=TreePos->Child_One;
     };
     if (TreePos->Type == MTK_COMMA) {
       if ( !TreePos->Parent || !TreePos->Child_One || !TreePos->Child_Two || (TreePos->Parent->Type != MTK_BESSJ))
         throw ERR_FUNCTION;
       temp = TreePos;
       TreePos = TreePos->Parent;
       TreePos->Child_One = temp->Child_One;
       TreePos->Child_Two = temp->Child_Two;
       temp->Child_One->Parent = TreePos; 
       temp->Child_Two->Parent = TreePos;
       delete temp;
     }
     else {
       if (TreePos->Child_Two==NULL) {
         TreePos->Brackets=1;
         TreePos=TreePos->Parent;
       }
       else {
         if (TreePos->Child_Two->Brackets == 1) {
           TreePos->Brackets=1;
           TreePos=TreePos->Parent;
         }
         else {
           TreePos=TreePos->Child_Two;
         }
       }
     }
   }
}


// Tidy up equation tree
void 
Parser::TidyUp()
{
  Node *TreePos;
  TreePos = equation->Root;
  
  equation->Reset(0);
  while ( equation->Root->Brackets == 0 ) {
    while ( TreePos->Child_One != NULL ) {
      if ( TreePos->Child_One->Brackets == 1 ) break;
      TreePos=TreePos->Child_One;
    };
    TreePos->Brackets = 1;
    TreePos = TreePos->Parent;
    if ( TreePos == NULL ) break;
    if ( TreePos->Child_Two == NULL ) {
      TidyNode(TreePos);
    }
    else {
      if ( (TreePos->Child_Two)->Brackets==0 ) {
        TreePos=TreePos->Child_Two;
      }
      else {
        TidyNode( TreePos );
      }
    }
  }
}


// Tidy up this node
void 
Parser::TidyNode(Node *&TreePos)
{
  Node *Temp;
  double value;

  switch (TreePos->Type) {
  case MTK_ADD: case MTK_SUBTRACT:
    if ((TreePos->Child_One->Type == MTK_NUMBER)&&(TreePos->Child_Two->Type == MTK_NUMBER)) {
      if (TreePos->Type == MTK_ADD) {
        TreePos->Value.Number = TreePos->Child_One->Value.Number + TreePos->Child_Two->Value.Number;
      }
      else {
        TreePos->Value.Number = TreePos->Child_One->Value.Number - TreePos->Child_Two->Value.Number;
      }
      TreePos->Type = MTK_NUMBER;
      TreePos->evalfn=func_details[MTK_NUMBER].func_ptr;
      equation->Delete(TreePos->Child_One);
      equation->Delete(TreePos->Child_Two);
      TreePos->Child_One = NULL;
      TreePos->Child_Two = NULL;
      break;
    }
    if ((TreePos->Child_One->Type == MTK_NUMBER) && ((TreePos->Child_One->Value).Number == 0)) {
      TreePos = equation->DeleteNode(TreePos,TreePos->Child_One,TreePos->Child_Two);
      break;
    }
    if ((TreePos->Child_Two->Type == MTK_NUMBER) && ((TreePos->Child_Two->Value).Number == 0)) {
      TreePos = equation->DeleteNode(TreePos,TreePos->Child_Two,TreePos->Child_One);
      break;
    };
    if ((TreePos->Child_Two->Type == MTK_NUMBER) && ((TreePos->Child_Two->Value).Number < 0)) {
      TreePos->Child_Two->Value.Number = -(TreePos->Child_Two->Value.Number);
      if (TreePos->Type == MTK_ADD) TreePos->Type = MTK_SUBTRACT;
      else TreePos->Type = MTK_ADD;
      TreePos->evalfn=func_details[TreePos->Type].func_ptr;
      break;
    };
    break;
  case MTK_MULTIPLY:
    if ((TreePos->Child_One->Type == MTK_NUMBER) && (TreePos->Child_Two->Type == MTK_NUMBER)) {
      TreePos->Value.Number = TreePos->Child_One->Value.Number * TreePos->Child_Two->Value.Number;
      TreePos->Type = MTK_NUMBER;
      TreePos->evalfn = func_details[MTK_NUMBER].func_ptr;
      equation->Delete(TreePos->Child_One);
      equation->Delete(TreePos->Child_Two);
      TreePos->Child_One = NULL;
      TreePos->Child_Two = NULL;
      break;
    }
    if (((TreePos->Child_One->Type == MTK_NUMBER) && ((TreePos->Child_One->Value).Number == 0)) 
        || ((TreePos->Child_Two->Type == MTK_NUMBER) && ((TreePos->Child_Two->Value).Number == 0))) {
      equation->Delete(TreePos->Child_One);
      equation->Delete(TreePos->Child_Two);
      TreePos->Type=MTK_NUMBER;
      TreePos->evalfn = func_details[MTK_NUMBER].func_ptr;
      (TreePos->Value).Number=0;
      TreePos->Child_One=NULL;
      TreePos->Child_Two=NULL;
      break;
    };
    if ((TreePos->Child_One->Type == MTK_NUMBER) && ((TreePos->Child_One->Value).Number == 1)) {
      TreePos = equation->DeleteNode(TreePos,TreePos->Child_One,TreePos->Child_Two);
      break;
    };
    if ((TreePos->Child_Two->Type == MTK_NUMBER) && ((TreePos->Child_Two->Value).Number == 1)) {
      TreePos = equation->DeleteNode(TreePos,TreePos->Child_Two,TreePos->Child_One);
      break;
    };
    if ((TreePos->Child_One->Type >= MTK_LN) && (TreePos->Child_One->Type <= MTK_FRAC)) {
      if ((TreePos->Child_Two->Type < MTK_LN) || (TreePos->Child_Two->Type > MTK_FRAC)) {
        Temp=TreePos->Child_One;
        TreePos->Child_One=TreePos->Child_Two;
        TreePos->Child_Two=Temp;
      };
    };
    break;
  case MTK_DIVIDE:
    if ((TreePos->Child_Two->Type == MTK_NUMBER) && ((TreePos->Child_Two->Value).Number == 0)) {
      throw ERR_DIV_ZERO;
    };
    if ((TreePos->Child_One->Type == MTK_NUMBER) && (TreePos->Child_Two->Type == MTK_NUMBER)) {
      TreePos->Value.Number = TreePos->Child_One->Value.Number / TreePos->Child_Two->Value.Number;
      TreePos->Type = MTK_NUMBER;
      TreePos->evalfn = func_details[MTK_NUMBER].func_ptr;
      equation->Delete(TreePos->Child_One);
      equation->Delete(TreePos->Child_Two);
      TreePos->Child_One = NULL;
      TreePos->Child_Two = NULL;
      break;
    }
    if ((TreePos->Child_One->Type == MTK_NUMBER) && ((TreePos->Child_One->Value).Number == 0)) {
      equation->Delete(TreePos->Child_One);
      equation->Delete(TreePos->Child_Two);
      TreePos->Type=MTK_NUMBER;
      TreePos->evalfn = func_details[MTK_NUMBER].func_ptr;
      (TreePos->Value).Number=0;
      TreePos->Child_One=NULL;
      TreePos->Child_Two=NULL;
      break;
    }
    if ((TreePos->Child_Two->Type == MTK_NUMBER) && ((TreePos->Child_Two->Value).Number == 1)) {
      TreePos = equation->DeleteNode(TreePos,TreePos->Child_Two,TreePos->Child_One);
      break;
    }
    break;
  case MTK_POWER :
    if (((TreePos->Child_One)->Type==MTK_NUMBER)&&(((TreePos->Child_One)->Value).Number==0)
        &&((TreePos->Child_Two)->Type==MTK_NUMBER)&&(((TreePos->Child_Two)->Value).Number==0)) {
      throw ERR_DIV_ZERO;
    };
    if (((TreePos->Child_One)->Type==MTK_NUMBER)&&(((TreePos->Child_One)->Value).Number==0)) {
      equation->Delete(TreePos->Child_One);
      equation->Delete(TreePos->Child_Two);
      TreePos->Child_One=NULL;
      TreePos->Child_Two=NULL;
      TreePos->Type=MTK_NUMBER;
      TreePos->evalfn = func_details[MTK_NUMBER].func_ptr;
      (TreePos->Value).Number=0;
      break;
    }
    if (((TreePos->Child_Two)->Type==MTK_NUMBER)&&(((TreePos->Child_Two)->Value).Number==0)) {
      equation->Delete(TreePos->Child_One);
      equation->Delete(TreePos->Child_Two);
      TreePos->Child_One=NULL;
      TreePos->Child_Two=NULL;
      TreePos->Type=MTK_NUMBER;
      TreePos->evalfn = func_details[MTK_NUMBER].func_ptr;
      (TreePos->Value).Number=1;
      break;
    }
    if (((TreePos->Child_Two)->Type==MTK_NUMBER)&&(((TreePos->Child_Two)->Value).Number==1)) {
      TreePos = equation->DeleteNode(TreePos, TreePos->Child_Two, TreePos->Child_One);
      break;
    }
    if (((TreePos->Child_One)->Type==MTK_NUMBER)&&((TreePos->Child_Two)->Type==MTK_NUMBER)) {
      TreePos->Value.Number = pow(TreePos->Child_One->Value.Number,TreePos->Child_Two->Value.Number);
      TreePos->Type = MTK_NUMBER;
      TreePos->evalfn = func_details[MTK_NUMBER].func_ptr;
      equation->Delete(TreePos->Child_One);
      equation->Delete(TreePos->Child_Two);
      TreePos->Child_One = NULL;
      TreePos->Child_Two = NULL;
      break;
    }
    break;
  case MTK_BESSJ:
  case MTK_COMBINATION:
  case MTK_PERMUTATION:
    if ((TreePos->Child_One->Type == MTK_NUMBER) && (TreePos->Child_Two->Type == MTK_NUMBER)) {
      switch (TreePos->Type) {
      case MTK_BESSJ:
        value = ansi_jn((int)TreePos->Child_One->Value.Number,TreePos->Child_Two->Value.Number);
        break;
      case MTK_COMBINATION:
        value = combination(TreePos->Child_One->Value.Number, TreePos->Child_Two->Value.Number);
        break;
      case MTK_PERMUTATION:
        value = permutation(TreePos->Child_One->Value.Number, TreePos->Child_Two->Value.Number);
        break;
      }
      TreePos->Type = MTK_NUMBER;
      TreePos->evalfn = func_details[MTK_NUMBER].func_ptr;
      TreePos->Value.Number = value;
      equation->Delete(TreePos->Child_One);
      TreePos->Child_One = NULL;
      equation->Delete(TreePos->Child_Two);
      TreePos->Child_Two = NULL;
    }
    break;
  default :
    if ( ( (TreePos->Type) >= MTK_LN ) && ( (TreePos->Type) <= MTK_FACTORIAL) ) {
      if (TreePos->Child_One->Type == MTK_NUMBER) {
        if ( equation->Radians || (TreePos->Type < MTK_SIN) || ( TreePos->Type > MTK_ARCTAN) ) {
          value = (*(TreePos->evalfn))(TreePos->Child_One->Value.Number);
        }
        else {
          if ( TreePos->Type <= MTK_COT )
            value = (*(TreePos->evalfn))(TreePos->Child_One->Value.Number*DEGTORAD);
          else
            value = (*(TreePos->evalfn))(TreePos->Child_One->Value.Number)/DEGTORAD;
          }
        TreePos->Type = MTK_NUMBER;
        TreePos->evalfn = func_details[MTK_NUMBER].func_ptr;
        TreePos->Value.Number = value;
        equation->Delete(TreePos->Child_One);
        TreePos->Child_One = NULL;
      }
    }
  };

  TreePos->Brackets=1;
}

// Convert tree to text
void
Parser::TreeToText()
{
  size_t Length;
  char *Buffer;

  equation->Reset(0);
  Buffer = TreeToText(equation->Root);

  Length=strlen(Buffer);
  if (EquationText)
    delete [] EquationText;
  EquationText=new char[Length+1];
  strcpy(EquationText,Buffer);
  delete [] Buffer;
}

// Convert tree from node to text
char *
Parser::TreeToText(Node *TreeRoot)
{
  bool Flag=true;
  int length;
  int extra;
  Node *TreePos;
  char *Buffer, *ch1, *ch2, *tempbuffer;
  Buffer = new char[1024];
  int BufferSize = 1023;
  TreePos = TreeRoot;
  Buffer[0] = 0;

  while (Flag) {
    while (TreePos->Child_One != NULL){
      if (TreePos->Child_One->Brackets == 1) break;
      if ((TreePos->Type >= MTK_LN) && (TreePos->Type <= MTK_FRAC)) {
        NodeText(TreePos->Type,TreePos,Buffer, BufferSize);
      }
      if (TreePos->Type == MTK_BESSJ) {
        if (!TreePos->Child_One || !TreePos->Child_Two) {
          throw ERR_FUNCTION;
        }
        ch1 = TreeToText(TreePos->Child_One);
        ch2 = TreeToText(TreePos->Child_Two);
        length = int(strlen(ch1) + strlen(ch2) + 9);
        extra = BufferSize - int(strlen(Buffer));
        if (extra < length) {
          tempbuffer = Buffer;
          Buffer = new char[BufferSize + extra + 100];
          strcpy(Buffer, tempbuffer);
          delete [] tempbuffer;
        }
        strcat(Buffer,"( ");
        strcat(Buffer,ch1);
        strcat(Buffer," , ");
        strcat(Buffer,ch2);
        strcat(Buffer," ) ");
        delete [] ch1;
        delete [] ch2;
        TreePos->Brackets = 1;
        break;
      }
      else
        if ( Priority(TreePos->Child_One->Type, TreePos->Type) == false ) {
          NodeText(MTK_LBRACKET, NULL, Buffer, BufferSize);
        };
      TreePos=TreePos->Child_One;
    };
    if ((TreePos->Type<MTK_LN)||(TreePos->Type>MTK_FRAC)) {
      NodeText(TreePos->Type,TreePos,Buffer, BufferSize);
    }
    if ((TreePos->Child_Two == NULL) || (TreePos->Type == MTK_BESSJ)) {
      do {
        if (TreePos != TreeRoot) {
          if ((Priority(TreePos->Type, TreePos->Parent->Type) == false)
              && (TreePos->Child_One != NULL)){
            NodeText(MTK_RBRACKET, NULL, Buffer, BufferSize);
          }
          if (TreePos->Parent->Child_Two == TreePos) {
            if (((TreePos->Parent->Type == MTK_SUBTRACT)
                 && ((TreePos->Type == MTK_ADD) || (TreePos->Type == MTK_SUBTRACT)))
                || ((TreePos->Parent->Type == MTK_DIVIDE)
                    && ((TreePos->Type == MTK_MULTIPLY) || (TreePos->Type == MTK_DIVIDE)))) {
              NodeText(MTK_RBRACKET, NULL, Buffer, BufferSize);
            }
          }
        }
        if (TreePos == TreeRoot) {
          Flag = false;
          break;
        };
        TreePos=TreePos->Parent;
      }
      while (TreePos->Brackets==1);
    }
    else {
      if ((Priority(TreePos->Child_Two->Type,TreePos->Type)==false) ||
          ((TreePos->Type == MTK_SUBTRACT) && ((TreePos->Child_Two->Type == MTK_ADD) || (TreePos->Child_Two->Type == MTK_SUBTRACT))) ||
          ((TreePos->Type == MTK_DIVIDE) && ((TreePos->Child_Two->Type == MTK_MULTIPLY) || (TreePos->Child_Two->Type == MTK_DIVIDE))))
        NodeText(MTK_LBRACKET, NULL, Buffer, BufferSize);
      TreePos = TreePos->Child_Two;
    }
  }
  return Buffer;
}

// Convert node to text
void
Parser::NodeText(int type, Node *TreePos,char *Buffer,int &BufferSize)
{
   char Temp[30];
   char *tempbuffer;
   int extra;
   if ((type <1) || (type>MTK_NUMBER_OF_FUNCS)) {
     if (TreePos) TreePos->Brackets=1;
     throw ERR_FUNCTION;
   }
   if ( (type != MTK_NUMBER) && (type != MTK_VARIABLE) ) {
     if ( (type != MTK_MULTIPLY) || ( ((TreePos->Child_One)->Type!=MTK_NUMBER) || ((TreePos->Child_Two)->Type==MTK_NUMBER) ) ) {
        extra = int(strlen(Buffer)) + int(strlen(func_details[type].func_text)) - BufferSize;
        if ( extra >= -1 ) {
          extra +=100;
          tempbuffer = Buffer;
          Buffer = new char[BufferSize + extra];
          strcpy(Buffer,tempbuffer);
          delete [] tempbuffer;
        }
        strcat(Buffer,func_details[type].func_text);
      }
   }
   else {
     if (type == MTK_NUMBER){
       sprintf(Temp,"%.15G",(TreePos->Value).Number);
     }
     else {
       sprintf(Temp,"%s",equation->Variables.at(TreePos->Value.WhichVar - 1)->Name);
     }
     extra = int(strlen(Buffer)) + int(strlen(Temp)) - BufferSize;
     if ( extra >= -1 ) {
       extra +=50;
       tempbuffer = Buffer;
       Buffer = new char[BufferSize + extra];
       strcpy(Buffer,tempbuffer);
       delete [] tempbuffer;
     }
     strcat(Buffer,Temp);
   }

   if (TreePos) TreePos->Brackets=1;
}









