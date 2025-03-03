#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main (int argc, char *argv[])
{
  long i,j,k,nnp,nns,nep,nes,nnep,nnes,sts,stp;
  long *anp,**nods,**nodp;
  double z,xo,yo;
  double *xs,*ys,*xp,*yp;
  FILE *in,*ins,*inp,*out;
  
  ins = fopen(argv[1],"rt");
  inp = fopen(argv[2],"rt");
  in  = fopen(argv[3],"rt");
  out = fopen(argv[4],"wt");
  
  // ********
  //  SOIL
  // ********
  
  //  number of nodes on soil
  fscanf (ins,"%ld",&nns);
  
  xs = new double [nns];
  ys = new double [nns];
  
  for (i=0;i<nns;i++){
    fscanf (ins,"%le %le %le %ld",xs+i,ys+i,&z,&j);
  }
  
  //  number of elements
  fscanf (ins,"%ld %ld",&nes,&sts);
  
  nods = new long* [nes];
  for (i=0;i<nes;i++){
    fscanf (ins,"%ld",&nnes);
    nods[i] = new long [nnes];
    for (j=0;j<nnes;j++){
      fscanf (ins,"%ld",&nods[i][j]);
    }
    fscanf (ins,"%ld",&j);
  }
  

  // ********
  //  PLATE
  // ********

  //  number of nodes on plate
  fscanf (inp,"%ld",&nnp);
  
  xp = new double [nnp];
  yp = new double [nnp];
  
  for (i=0;i<nnp;i++){
    fscanf (inp,"%le %le %le %ld",xp+i,yp+i,&z,&j);
  }
  
  //  number of elements
  fscanf (inp,"%ld %ld",&nep,&stp);
  
  if (sts!=stp){
    fprintf (stderr,"\n\n different types of elements in soil and plate files \n");
  }

  nodp = new long* [nep];
  for (i=0;i<nep;i++){
    fscanf (inp,"%ld",&nnep);
    nodp[i] = new long [nnep];
    for (j=0;j<nnep;j++){
      fscanf (inp,"%ld",&nodp[i][j]);
    }
    fscanf (inp,"%ld",&j);
  }
  
  
  //  coordinates of origin of plate
  fscanf (in,"%le %le",&xo,&yo);
  
  
  //
  //  recalculation of coordinates of plate
  //
  for (i=0;i<nnp;i++){
    xp[i]+=xo;
    yp[i]+=yo;
  }
  

  //  auxiliary array of plate node numbers
  anp = new long [nnp];
  
  //  plate node correction
  for (i=0;i<nnp;i++){
    for (j=0;j<nns;j++){
      if (xp[i]==xs[j] && yp[i]==ys[j]){
	anp[i]=j+1;
	break;
      }
    }
  }
  
  //  
  for (i=0;i<nep;i++){
    for (j=0;j<nnep;j++){
      k=nodp[i][j]-1;
      nodp[i][j]=anp[k];
    }
  }
  
  
  //
  //  print of one file with mesh
  //
  
  //  nodes
  fprintf (out,"%ld\n",nns);
  for (i=0;i<nns;i++){
    fprintf (out,"%le %le 0.0 0\n",xs[i],ys[i]);
  }
  
  //  elements
  fprintf (out,"%ld %ld\n",nes+nep,sts);
  
  //  soil
  for (i=0;i<nes;i++){
    fprintf (out,"%ld  ",nnes);
    for (j=0;j<nnes;j++){
      fprintf (out,"  %ld",nods[i][j]);
    }
    fprintf (out,"  0\n");
  }
  
  //  plate
  for (i=0;i<nep;i++){
    fprintf (out,"%ld  ",nnep);
    for (j=0;j<nnep;j++){
      fprintf (out,"  %ld",nodp[i][j]);
    }
    fprintf (out,"  1\n");
  }
  
  
}
