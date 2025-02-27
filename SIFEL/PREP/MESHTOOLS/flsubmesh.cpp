#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*
  code generates a file with mesh of reinforcement-matrix problem
  the code merges several files with meshes together
  
*/
int main (int argc, char *argv[])
{
  long i,j,k,l,m;
  long nf,nn,ne,nln;
  long *cnns,*nns,*nes,**nne,**np,**ep,***nodes,*matmultip,*fibnum,*fibnod,*aux;
  double xx,yy,zz;
  double **x,**y,**z;
  FILE *in,*out;
  char inname[BUFSIZ],infile[BUFSIZ];
  char suffixtop[]=".top";

  if (argc!=4){
    printf("Wrong number of parameters!\n");
    printf ("flsubmesh <mesh file name> <number of fibres> <output>\n");
    exit(EXIT_FAILURE);
  }
  
  //  inname contains name of the files with meshes
  //  it has to be without suffix, e.g. top, etc., and without number
  strcpy(inname,argv[1]);
  
  //  number of fibres
  nf=atoi(argv[2]);

  
  //  array containing number of nodes on subdomains
  nns = new long [nf+1];
  //  array containing number of elements on subdomains
  nes = new long [nf+1];
  //  array containing cumulative numbers of nodes on subdomains
  cnns = new long [nf+1];
  
  //  array for nodal coordinates
  x = new double* [nf+1];
  y = new double* [nf+1];
  z = new double* [nf+1];

  //  nodal properties
  np = new long* [nf+1];
  //  element properties
  ep = new long* [nf+1];
  //  number of nodes on element
  nne = new long* [nf+1];
  //  nodes on elements
  nodes = new long** [nf+1];
  
  // *************************
  //  loop over input files
  //  reading of input files
  // *************************
  for (i=0;i<nf+1;i++){

    //  file with mesh of the i-th subdomain
    sprintf(infile,"%s%d",inname,i);
    strcat(infile,suffixtop);	
    in=fopen(infile,"r");
    if (in==NULL){
      fprintf (stderr,"\n File with the mesh of the i-th subdomain cannot be opened\n");
      exit(EXIT_FAILURE);
    }
    
    //  number of nodes on the i-th subdomain
    fscanf (in,"%ld",nns+i);
    
    //  node coordinates
    x[i] = new double [nns[i]];
    y[i] = new double [nns[i]];
    z[i] = new double [nns[i]];
    //  nodal properties
    np[i] = new long [nns[i]];
    
    for (j=0;j<nns[i];j++){
      fscanf (in,"%lf %lf %lf  %ld",&x[i][j],&y[i][j],&z[i][j],&np[i][j]);
    }
    
    //  number of elements on the i-th subdomain
    //  type of elements
    fscanf (in,"%ld %ld",nes+i,&j);
    
    //  element properties
    ep[i] = new long [nes[i]];
    //  nodes on elements
    nodes[i] = new long* [nes[i]];
    //  number of nodes on element
    nne[i] = new long [nes[i]];
    
    for (j=0;j<nes[i];j++){
      fscanf (in,"%ld",&nne[i][j]);
      nodes[i][j] = new long [nne[i][j]];
      for (k=0;k<nne[i][j];k++){
	fscanf (in,"%ld",&nodes[i][j][k]);
      }
      fscanf (in,"%ld",&ep[i][j]);
    }
    
  }
  fclose (in);
  
  
  // *********************************************************************
  //  determination of contact nodes between composite matrix and fibres
  // *********************************************************************
  //  array of multiplicity of composite matrix nodes
  matmultip = new long [nns[0]];
  //  array of fibre numbers
  //  fibnum[i]=j - the i-th node of composite matrix is connected with the j-th fibre
  fibnum = new long [nns[0]];
  //  array of fibre nodes connected to composite matrix nodes
  //  fibnod[i]=k - the i-th node of composite matrix is connected with the k-th node (of the j-th fibre)
  fibnod = new long [nns[0]];
  
  //  loop over nodes on composite matrix
  for (i=0;i<nns[0];i++){
    matmultip[i]=0;
    
    //  coordinates of actual node
    xx=x[0][i];
    yy=y[0][i];
    zz=z[0][i];
    
    //  loop over fibres
    for (j=1;j<nf+1;j++){
      for (k=0;k<nns[j];k++){
	if (xx==x[j][k] && yy==y[j][k] && zz==z[j][k]){
	  matmultip[i]++;
	  //  number of connected fibre
	  fibnum[i]=j;
	  //  number of connected node
	  fibnod[i]=k;
	}
      }
    }
    
    if (matmultip[i]!=0 && matmultip[i]!=1){
      fprintf (stderr,"\n the %ld-th node of composite matrix is connected to %ld fibre nodes\n",i,matmultip[i]);
      exit(EXIT_FAILURE);
    }
  }
  
  
  // **************************************************
  //  merging of meshes, node and element renumbering
  // **************************************************
  out = fopen (argv[3],"w");
  
  nn=0;  ne=0;
  for (i=0;i<nf+1;i++){
    nn+=nns[i];
    ne+=nes[i];
  }
  
  //  total number of nodes
  fprintf (out,"%ld\n",nn);
  
  for (i=0;i<nf+1;i++){
    for (j=0;j<nns[i];j++){
      fprintf (out,"%20.15le %20.15le %20.15le  %ld\n",x[i][j],y[i][j],z[i][j],np[i][j]);
    }
  }
  
  //  number of elements
  fprintf (out,"%ld %ld\n",ne,3);

  l=0;
  for (i=0;i<nf+1;i++){
    aux = new long [nns[i]];
    
    for (j=0;j<nns[i];j++){
      aux[j]=j+l+1;
    }
    
    for (j=0;j<nes[i];j++){
      //  number of nodes on element
      fprintf (out,"%ld  ",nne[i][j]);
      //  node numbers on element
      for (k=0;k<nne[i][j];k++){
	m=aux[nodes[i][j][k]-1];
	fprintf (out,"%ld  ",m);
      }
      //  element property
      fprintf (out,"%ld\n",ep[i][j]);
    }
    
    delete [] aux;
    l+=nns[i];
  }
  
  //  number of connected (layered) nodes
  nln=0;
  for (i=0;i<nns[0];i++){
    if (matmultip[i]==1){
      nln++;
    }
  }
  
  cnns[0]=0;
  for (i=1;i<nf+1;i++){
    cnns[i]=cnns[i-1]+nns[i-1];
  }
  
  fprintf (out,"%ld\n",nln);
  //  list of connected nodes
  j=0;
  for (i=0;i<nns[0];i++){
    if (matmultip[i]==1){
      k=cnns[fibnum[i]]+fibnod[i]+1;
      fprintf (out,"%ld  2  %ld %ld\n",j+1,i+1,k);
      j++;
    }
  }

  fprintf (out,"\n");
  fclose (out);
}
