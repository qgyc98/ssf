/*
  File:     qen2delem.cpp
  Author:   Tomas Koudelka, 3.2011
  Purpose:  Fictious general 2D element for HERMES adaptivity
*/

#include "globalt.h"
#include "gen2delem.h"
#include "genfile.h"
#include "globmatt.h"

gen2delem::gen2delem (void)
{
  long i;
  
  //  number of nodes on element
  nne=1;
  //  number of edges
  ned=0;
  //  number of nodes on one edge
  nned=0;
  //  geometrical dimension
  ncomp=2;
  
  //  number of transported variables
  ntm=Tp->ntm;
  
  
  nip = new long* [ntm];
  for (i=0;i<ntm;i++){
    nip[i] = new long [ntm];
  }
  
  
  switch (Tp->tmatt){
  case nomedium:{  break; }
  case onemedium:{
    nip[0][0]=1;
    ndofe=1;
    break;
  }
  case twomediacoup:{
    if (Tp->savemode==0){
      nip[0][0]=1;       nip[0][1]=1;       nip[1][0]=1;       nip[1][1]=1;
    }
    if (Tp->savemode==1){
      nip[0][0]=1;       nip[0][1]=0;       nip[1][0]=0;       nip[1][1]=0;
    }
    ndofe=2;
    break;
  }
  case threemediacoup:{
    if (Tp->savemode==0){
      nip[0][0]=1;  nip[0][1]=1;  nip[0][2]=1;
      nip[1][0]=1;  nip[1][1]=1;  nip[1][2]=1;
      nip[2][0]=1;  nip[2][1]=1;  nip[2][2]=1;
    }
    if (Tp->savemode==1){
      nip[0][0]=1;  nip[0][1]=0;  nip[0][2]=0;
      nip[1][0]=0;  nip[1][1]=0;  nip[1][2]=0;
      nip[2][0]=0;  nip[2][1]=0;  nip[2][2]=0;
    }
    ndofe=3;
    break;
  }
  case fourmediacoup:{
    if (Tp->savemode==0){
      nip[0][0]=1;  nip[0][1]=1;  nip[0][2]=1;  nip[0][3]=1;
      nip[1][0]=1;  nip[1][1]=1;  nip[1][2]=1;  nip[1][3]=1;
      nip[2][0]=1;  nip[2][1]=1;  nip[2][2]=1;  nip[2][3]=1;
      nip[3][0]=1;  nip[3][1]=1;  nip[3][2]=1;  nip[3][3]=1;
    }
    if (Tp->savemode==1){
      nip[0][0]=1;  nip[0][1]=0;  nip[0][2]=0;  nip[0][3]=0;
      nip[1][0]=0;  nip[1][1]=0;  nip[1][2]=0;  nip[1][3]=0;
      nip[2][0]=0;  nip[2][1]=0;  nip[2][2]=0;  nip[2][3]=0;
      nip[3][0]=0;  nip[3][1]=0;  nip[3][2]=0;  nip[3][3]=0;
    }
    ndofe=4;
    break;
  }
  default:{
    print_err("unknown number of transported matters is required",__FILE__,__LINE__,__func__);
  }
  }
}

gen2delem::~gen2delem (void)
{
  long i;

  for (i=0;i<ntm;i++){
    delete [] nip[i];
  }
  delete [] nip;
}

