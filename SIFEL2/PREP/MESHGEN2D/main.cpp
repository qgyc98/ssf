#include <stdio.h>
#include <stdlib.h>
#include "siftop.h"
#include "meshgen2d.h"

siftop Top;

int main (int argc,char *argv[])
{
  XFILE *inedges;
  XFILE *in;
  FILE *out;
  long j=0;
  long i,w,z,l,k=0;
  long **edges;
  long ***edgenod;
  long *edgdiv;       // pole, kde je ulozeno deleni kazde krajni hrany v jedne topologii(rezu)
  long edgcount;      // pocet hran v makroelementu
  long **elemdiv;
  siftop rez1;


  fprintf(stdout, "\n\n\n---          GENERATOR OF STRUCTURED MESH             ---\n");
  if (argc < 4){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : %s input_file_name input_edges_file_name output_file_name \n\n", argv[0]);
    return(1);
  }
  out = fopen(argv[3], "wt");
  if (out==NULL)
  {
    fprintf (stderr,"\n Output file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }
  
  
  in = xfopen(argv[1], "rt");
  if (in==NULL)
  {
    fprintf (stderr,"\n Input file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }

  inedges = xfopen(argv[2], "rt");
  if (in==NULL)
  {
    fprintf (stderr,"\n Input file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }

   read_edgdiv(inedges,edgdiv,edgcount);
   Top.read(in, 0, 1);

  edges = new long*[4];
  for (i=0;i<4;i++)
  {
    edges[i]= new long[3];
    for (z=0;z<3;z++) edges [i][z]=0;
  }

  edgenod = new long**[Top.ne];
  for (i=0;i<Top.ne;i++)
  {
    edgenod[i]= new long*[4];
    for (w=0;w<4;w++)
    {
      edgenod[i][w] = new long [11];
      for (z=0;z<11;z++) edgenod[i][w][z]=-1;
    }
  }

  elemdiv = new long*[Top.ne];
  for (i=0;i<Top.ne;i++)
  {
    elemdiv[i]= new long[4];
    for (z=0;z<4;z++) elemdiv [i][z]=0;
  }

  rez1.ne=0;rez1.nn=0;  
  for (i=0;i<Top.ne;i++)
  {
    searchsimedg(Top.elements,i,edges);
    div_elements (edgdiv,elemdiv,Top.elements[i].propedg,edges,i);
    rez1.ne += elemdiv[i][0]*elemdiv[i][1];
    rez1.nn += (elemdiv[i][0]+1)*(elemdiv[i][1]+1);
  }
  
  rez1.nodes = new snode[rez1.nn];
  rez1.elements = new selement[rez1.ne];
  for (i=0;i<rez1.ne;i++)
  {
    rez1.elements[i].type = isolinear2d;
    rez1.elements[i].alloc(1);
  }

  for (i=0;i<Top.ne;i++)
  {
    l=j;
    searchsimedg(Top.elements,i,edges);
    gennodes(elemdiv[i][0],elemdiv[i][1],Top.nodes,Top.elements[i].nodes,rez1.nodes,j,edges,edgenod,i,Top.elements[i]);
    genelements(elemdiv[i][0],elemdiv[i][1],Top.elements[i],rez1.elements,edgenod[i],edges,l,k);
  }

 
  rez1.print(out);
  fclose(out);

  fprintf(stdout, "\nMesh generation has been finished\n");
  fprintf(stdout, "Number of nodes: %ld\n""Number of elements: %ld\n", rez1.nn, rez1.ne);
  return 0;
}

// nacist deleni zvlast ze souboru
