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

#include "tree.h"
#include "functext.h"

TTree::TTree()
{
   Root=new Node;
   Root->Parent=NULL;
   Root->Child_One=NULL;
   Root->Child_Two=NULL;
   Root->Brackets=0;
   Root->evalfn = NULL;
}

TTree::~TTree()
{
   Delete(Root);
}

void
TTree::InitNode(Node *Initnode, func_type type, NumVar val, int brack)
{
  Initnode->Type=type;
  Initnode->Value=val;
  Initnode->Brackets=brack;
  Initnode->evalfn = func_details[type].func_ptr;
}

void
TTree::InitNode(Node *Initnode, func_type type)
{
  Initnode->Type = type;
  Initnode->Value.Number = 0.0;
  Initnode->Brackets = 0;
  Initnode->evalfn = func_details[type].func_ptr;
}

Node *
TTree::AddParent(Node *TreePos)
{
   Node *Temp;
   Temp=TreePos->Parent;
   (TreePos->Parent)=new Node;
   (TreePos->Parent)->Child_One=TreePos;
   TreePos=TreePos->Parent;
   TreePos->Parent=Temp;
   if (Temp!=NULL)
   {
      if (Temp->Child_One==TreePos->Child_One)
      {
         Temp->Child_One=TreePos;
      }
      else
      {
         Temp->Child_Two=TreePos;
      };
   }
   else
   {
      Root=TreePos;
   };
   return TreePos;
}

Node *
TTree::AddChildOne(Node *TreePos)
{
   Node *Temp;
   Temp=TreePos->Child_One;
   TreePos->Child_One=new Node;
   (TreePos->Child_One)->Parent=TreePos;
   TreePos=TreePos->Child_One;
   TreePos->Child_One=Temp;
   if (Temp!=NULL)
   {
      Temp->Parent=TreePos;
   };
   return TreePos;
}

Node *
TTree::AddChildTwo(Node *TreePos)
{
   Node *Temp;
   Temp=TreePos->Child_Two;
   TreePos->Child_Two=new Node;
   (TreePos->Child_Two)->Parent=TreePos;
   TreePos=TreePos->Child_Two;
   TreePos->Child_One=Temp;
   if (Temp!=NULL) Temp->Parent=TreePos;
   return TreePos;
}

Node *
TTree::DeleteNode(Node *TreePos, Node *Blank, Node *Other)
{
  Node *Temp;
  if (Blank==NULL) {
    Blank=TreePos->Child_Two;
    Other=TreePos->Child_One;
  };
  if (Blank!=NULL) Delete(Blank);
  if (TreePos==Root) {
    Root=Other;
    Root->Parent=NULL;
    delete TreePos;
    return Root;
  };
  if (Other!=NULL) {
    Other->Parent=TreePos->Parent;
  };
  Temp=TreePos;
  TreePos=TreePos->Parent;
  if (TreePos!=NULL) {
    if (TreePos->Child_One==Temp) TreePos->Child_One=Other;
    else TreePos->Child_Two=Other;
  }
  delete Temp;
  return Other;
}

void 
TTree::Reset(int SetTo,Node *Start)
{
   if (Start==NULL) Start=Root;
   Node *TreePos;
   TreePos=Start;
   while (Start->Brackets!=SetTo)
   {
       while (TreePos->Child_One!=NULL)
       {
	  if ((TreePos->Child_One)->Brackets==SetTo) break;
          TreePos=TreePos->Child_One;
       };
       TreePos->Brackets=SetTo;
       if (TreePos==Start) break;
       TreePos=TreePos->Parent;
       if (TreePos->Child_Two!=NULL)
       {
	  if ((TreePos->Child_Two)->Brackets!=SetTo)
	  {
	     TreePos=TreePos->Child_Two;
	  }
	  else
	  {
	     TreePos->Brackets=SetTo;
          };
       }
       else
       {
	  TreePos->Brackets=SetTo;
       };
   };
}

Node *
TTree::Copy(Node *From,Node *To)
{
   Node *PointerFrom=From;
   Node *PointerTo=To;

   if (PointerFrom==NULL) return NULL;
   To->Brackets=0;
   PointerTo->Child_Two=NULL;
   PointerTo->Type=PointerFrom->Type;
   PointerTo->evalfn = PointerFrom->evalfn;
   PointerTo->Value=PointerFrom->Value;
   PointerTo->CurrentValue=0;
   while (To->Brackets==0)
   {
      while (PointerFrom->Child_One!=NULL)
      {
	 PointerTo->Child_One=new Node;
         (PointerTo->Child_One)->Parent=PointerTo;
         PointerTo=PointerTo->Child_One;
	 PointerFrom=PointerFrom->Child_One;
         PointerTo->Child_Two=NULL;
         PointerTo->Type=PointerFrom->Type;
	 PointerTo->evalfn=PointerFrom->evalfn;
	 PointerTo->Value=PointerFrom->Value;
         PointerTo->Brackets=0;
	 PointerTo->CurrentValue=0;
      };
      PointerTo->Child_One=NULL;
      PointerTo->Brackets=1;
      while (PointerFrom!=From)
      {
	 if ((PointerFrom->Child_Two!=NULL)&&(PointerTo->Child_Two==NULL)) break;
         PointerFrom=PointerFrom->Parent;
         PointerTo->Brackets=1;
         PointerTo=PointerTo->Parent;
      };
      if ((PointerFrom==From)&&((PointerFrom->Child_Two==NULL)
	||(PointerTo->Child_Two!=NULL))) break;
      if ((PointerFrom->Child_Two!=NULL)&&(PointerTo->Brackets==0))
      {
         PointerTo->Child_Two=new Node;
         (PointerTo->Child_Two)->Parent=PointerTo;
         PointerTo=PointerTo->Child_Two;
	 PointerFrom=PointerFrom->Child_Two;
         PointerTo->Brackets=0;
	 PointerTo->Type=PointerFrom->Type;
         PointerTo->evalfn=PointerFrom->evalfn;
         PointerTo->Value=PointerFrom->Value;
	 PointerTo->CurrentValue=0;
         PointerTo->Child_Two=NULL;
      }
      else
      {
	 To->Brackets=1;
      };
   };
   Reset(2,To);
   Reset(0,To);
   return To;
}

void 
TTree::Delete(Node *From)
{
  Node *Temp=From;
  Node *Marked=NULL;
  if (From==NULL) return;
  while ((From->Child_One!=NULL)||(From->Child_Two!=NULL)) {
    while (Temp->Child_One!=NULL) {
      Temp=Temp->Child_One;
    };
    Marked=Temp;
    Temp=Temp->Parent;
    if (Temp->Child_One==Marked) {
      Temp->Child_One=NULL;
    }
    else {
      Temp->Child_Two=NULL;
    };
    delete Marked;
    while ((Temp->Child_Two==NULL)&&(Temp!=From)) {
      Marked=Temp;
      Temp=Temp->Parent;
      if (Temp->Child_One==Marked) {
	Temp->Child_One=NULL;
      }
      else {
	Temp->Child_Two=NULL;
      };
      delete Marked;
      Temp->Child_One=NULL;
    };
    if (Temp->Child_Two!=NULL) {
      Temp=Temp->Child_Two;
    }
    else {
      break;
    };
  };
  delete From;
}


Node *
TTree::AddBelowOne(Node *Pointer,func_type Type,double Value)
{
   Pointer->Child_One=new Node;
   (Pointer->Child_One)->Parent=Pointer;
   Pointer=Pointer->Child_One;
   Pointer->Child_One=NULL;
   Pointer->Child_Two=NULL;
   Pointer->Type=Type;
   Pointer->evalfn = func_details[Type].func_ptr;
   (Pointer->Value).Number=Value;
   Pointer->Brackets=0;
   return Pointer;
}

Node *
TTree::AddBelowTwo(Node *Pointer,func_type Type,double Value)
{
   Pointer->Child_Two=new Node;
   (Pointer->Child_Two)->Parent=Pointer;
   Pointer=Pointer->Child_Two;
   Pointer->Child_One=NULL;
   Pointer->Child_Two=NULL;
   Pointer->Type=Type;
   Pointer->evalfn = func_details[Type].func_ptr;
   (Pointer->Value).Number=Value;
   Pointer->Brackets=0;
   return Pointer;
}
