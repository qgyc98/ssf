#include <stdio.h>
#include <stdlib.h>
#include "siftop.h"
#include "bridgen.h"

//#include "gtopology.h"
//#include "difcalc.h"

siftop Top; // Vysledna topologie

int main (int argc,char *argv[])
{
  XFILE *in;
  FILE *out;
  long j=0;
  long i,w,z,l,h,k=0;
  long **edges;
  long **corner;
  long ****edgenod;
  long *edgdiv;       // pole, kde je ulozeno deleni kazde krajni hrany v jedne topologii(rezu)
  long ***elemdiv;
  long numcuts;       // pocet char. rezu
  long numel=0,numnod=0;
  long numsec,dn;
  double zetko=0.0;
  double *pom;
  section *arsec;
  siftop *cuts;       // pole siftopu, kde jsou ulozeny zadane charakteristicke rezy
  siftop *topcuts;    // pole siftopu, do ktereho jsou vygenerovany charakteristicke rezy

 

  fprintf(stdout, "\n\n\n---          GENERATOR OF STRUCTURED MESH             ---\n");
  fprintf(stdout, "                  version 1.1  \n\nCreated by Josef Fiedler 09/2012\n\n\n");
  if (argc < 3){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : %s input_file_name output_file_name \n\n", argv[0]);
    return(1);
  }
  in = xfopen(argv[1], "rt");
  if (in==NULL)
  {
    fprintf (stderr,"\n input file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }
  out = fopen(argv[2], "wt");
  if (out==NULL)
  {
    fprintf (stderr,"\n Output file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }

  read (in,edgdiv,arsec,cuts,numcuts,topcuts,numsec); // cteni ze souboru
  xfclose(in);
  alloc (edges,elemdiv,numcuts,cuts,corner);  // alokovani edges a elemdiv


/*
Nasledujic cyklus:
- spocita, kolik bude uzlu a elementu v nasledne vygenerovanych charakteristickych rezech
- naplni pole elemdiv, ktere udava deleni na hranach v makroelementech z jednotlivych char. rezu
*/
  for (w=0;w<numcuts;w++)
  {
    for (i=0;i<cuts[w].ne;i++)
    {
      searchsimedg(cuts[w].elements,i,edges,corner);
      div_elements (edgdiv,elemdiv[w],cuts[w].elements[i].propedg,edges,i);
      topcuts[w].ne += elemdiv[w][i][0]*elemdiv[w][i][1];
      topcuts[w].nn += (elemdiv[w][i][0]+1)*(elemdiv[w][i][1]+1);
    }
  }

  edgenod = new long***[numcuts];
  for (h=0;h<numcuts;h++)
  {
    edgenod[h] = new long**[cuts[h].ne];
    for (i=0;i<cuts[h].ne;i++)
    {
      edgenod[h][i]= new long*[4];
      for (w=0;w<4;w++)
      {
        edgenod[h][i][w] = new long [elemdiv[h][i][w]+1];
        for (z=0;z<(elemdiv[h][i][w]+1);z++) edgenod[h][i][w][z]=-1;
      }
    }
  }
  for (h=0;h<numcuts;h++)               // alokovani pole siftopu topcuts
  {
    topcuts[h].nodes = new snode[topcuts[h].nn];
    topcuts[h].elements = new selement[topcuts[h].ne];
    for (i=0;i<topcuts[h].ne;i++)
    {
      topcuts[h].elements[i].type = isolinear2d;
      topcuts[h].elements[i].alloc(1);
    }
  }

  

  printf("\n\n");
  printf("Generovani charakteristickych rezu: \n");
  
  for (w=0;w<numcuts;w++)           // vygenerovani vsech char. rezu/naplneni pole topcuts
  {
    printf(" - Rez %ld ... ",w+1);
    l=0;j=0;k=0;
    for (i=0;i<cuts[w].ne;i++)
    {
      l=j;
      searchsimedg(cuts[w].elements,i,edges,corner);
      gennodes(elemdiv[w][i][0],elemdiv[w][i][1],cuts[w].nodes,cuts[w].elements[i].nodes,topcuts[w].nodes,j,edges,corner,edgenod[w],elemdiv[w],i);
      genelements(elemdiv[w][i][0],elemdiv[w][i][1],cuts[w].elements[i].prop,cuts[w].elements[i].propedg,cuts[w].elements[i].propsurf,topcuts[w].elements,edgenod[w][i],edges,corner,l,k);
    }
    topcuts[w].nn=j;
    printf("OK \n");
  }

  // check consistency of the beginning and end sections, i.e. number of nodes and elements must be same.
  long ret=0L;
  for (i=0;i<numsec;i++)
  {
    if ((topcuts[arsec[i].a].nn != topcuts[arsec[i].b].nn) ||
        (topcuts[arsec[i].a].ne != topcuts[arsec[i].b].ne))
    {
      fprintf(stderr, "\nSegment %ld has incompatible meshes at the beginning and end sections.\n"
              " - beginning section %ld: nn = %ld, ne = %ld.\n"
              " - end section %ld: nne = %ld, ne = %ld.\n", i+1,
              arsec[i].a+1, topcuts[arsec[i].a].nn, topcuts[arsec[i].a].ne,
              arsec[i].b+1, topcuts[arsec[i].b].nn, topcuts[arsec[i].b].ne);
      ret = 1L;
    }
  }
    
  if(ret){
    return 1;
  }

  Top.ne=0;Top.nn=0;
  Top.nn+=topcuts[arsec[0].a].nn;
  for (i=0;i<numsec;i++)
  {    
    Top.ne += (arsec[i].division)*topcuts[arsec[i].a].ne;
    Top.nn += (arsec[i].division)*topcuts[arsec[i].a].nn;
    if (i > 0){
      dn = topcuts[arsec[i].a].nn - topcuts[arsec[i-1].a].nn;
      if (dn > 0)
        Top.nn += dn;
    }
  }
  
  pom = new double[Top.nn];  
  Top.nodes = new snode[Top.nn];
  Top.elements = new selement[Top.ne];
  for (i=0;i<Top.ne;i++)
  {
    Top.elements[i].type = isolinear3d;
    Top.elements[i].alloc(1);
  }


  printf("\n\n");
  printf("Generovani useku: \n");

  // take all inital section nodes from the section arsec.a
  copy_section_nodes2top(topcuts[arsec[0].a], 0, topcuts[arsec[0].a].nn, numnod, Top);
  
  long pcut_id = -1;
  for (i=0;i<numsec;i++)     // generovani useku
  {
    printf(" - Usek %ld ... ",i+1);
    if (i > 0){
      dn = topcuts[arsec[i].a].nn - topcuts[arsec[i-1].a].nn;
      if (dn > 0){
        copy_shift_set_section_nodes2top(topcuts[arsec[i].a], topcuts[arsec[i-1].a].nn, dn, 0.0, 0.0, zetko, numnod, Top);
      }
    }
    gensection2(pcut_id, arsec[i], topcuts, Top, numel, numnod, zetko);
    pcut_id = arsec[i].a;
    printf("OK \n");
  }
  printf("\n\n");
  printf("Celkovy pocet uzlu: %ld \n",Top.nn);
  printf("Celkovy pocet prvku: %ld \n",Top.ne);

/*
  gtopology gt;
  vector xx(8), yy(8), zz(8);
  double v;

  Top.exporttop(&gt);
  for (i=0;i<Top.ne;i++) 
  {
    gt.give_node_coord3d (xx, yy, zz, i);
    jac_3d (v, xx, yy, zz, 1.0, 1.0, 1.0);
    if (v < 0.0)
      fprintf(stdout, "Element %ld has negative jacobian\n", i+1);
  }
*/

  for (i=0;i<Top.nn;i++)    // prohozeni souradnic: x ve smeru strednice, y a z v rovine prurezu
  {
    pom[i]         = Top.nodes[i].y;
    Top.nodes[i].y = Top.nodes[i].x;
    Top.nodes[i].x = Top.nodes[i].z;
    Top.nodes[i].z = pom[i];
  }
  // uprava poradi uzlu na prvku a jeho hranovych a plosnych property
  // s ohledem na zaporny jakobian
  long aux;
  for (i=0; i<Top.ne; i++){
    // swap nodes 1 and 4
    aux = Top.elements[i].nodes[0];
    Top.elements[i].nodes[0] = Top.elements[i].nodes[3];
    Top.elements[i].nodes[3] = aux;
    // swap nodes 2 and 3
    aux = Top.elements[i].nodes[1];
    Top.elements[i].nodes[1] = Top.elements[i].nodes[2];
    Top.elements[i].nodes[2] = aux;
    // swap nodes 5 and 8
    aux = Top.elements[i].nodes[4];
    Top.elements[i].nodes[4] = Top.elements[i].nodes[7];
    Top.elements[i].nodes[7] = aux;
    // swap nodes 6 and 7
    aux = Top.elements[i].nodes[5];
    Top.elements[i].nodes[5] = Top.elements[i].nodes[6];
    Top.elements[i].nodes[6] = aux;

    // swap surface properties on surfaces 2 and 4
    aux = Top.elements[i].propsurf[1];
    Top.elements[i].propsurf[1] = Top.elements[i].propsurf[3];
    Top.elements[i].propsurf[3] = aux;
    
    // swap edge properties on edges 1 and 3
    aux = Top.elements[i].propedg[0];
    Top.elements[i].propedg[0] = Top.elements[i].propedg[2];
    Top.elements[i].propedg[2] = aux;
    // swap edge properties on edges 6 and 7
    aux = Top.elements[i].propedg[5];
    Top.elements[i].propedg[5] = Top.elements[i].propedg[6];
    Top.elements[i].propedg[6] = aux;
    // swap edge properties on edges 5 and 8
    aux = Top.elements[i].propedg[4];
    Top.elements[i].propedg[4] = Top.elements[i].propedg[7];
    Top.elements[i].propedg[7] = aux;
    // swap edge properties on edges 9 and 11
    aux = Top.elements[i].propedg[8];
    Top.elements[i].propedg[8] = Top.elements[i].propedg[10];
    Top.elements[i].propedg[10] = aux;    
  }

  // adding nodal properties according to element properties
  selement *auxe;
  snode    *auxn;
  for(i=0; i<Top.ne; i++)
  {
    auxe = Top.elements+i;
    for(j=0; j<auxe->nne; j++)
    {
      auxn = Top.nodes + auxe->nodes[j];
      if (auxn->nprop == 0)
      {
        auxn->nprop = 1;
        auxn->entid = new entityp[1];
        auxn->prop  = new long[1];
        auxn->entid[0] = eregion;
        auxn->prop[0]  = auxe->prop;
      }
    }
  }
  for(i=0; i<Top.ne; i++)
  {
    auxe = Top.elements+i;
    // property ploch se mohou s ohledem na zpusob generace 3D prvku
    // vyskytovat pouze u ploch c. 1-4
    for(j=0; j<auxe->nsurf/2; j++)
    {
      if (auxe->propsurf[j] != 0)
      {
        for(k=0; k<auxe->nnsurf[j]; k++)
        {
          auxn = Top.nodes + auxe->nodes[auxe->surfnod[j][k]];
          auxn->add_prop(esurface, auxe->propsurf[j]);
        }
      }
    }
  }

  printf("\n\n");
  printf("Vystup topologie ...");
  
  Top.print(out); // vystup
  printf("OK \n");
  fclose(out);
  printf("\nPress Enter to continue ...\n");
  scanf("%*c");

  delete [] edgdiv;
  dealloc (edges,elemdiv,numcuts,cuts,corner);

  for (h=0;h<numcuts;h++){
    for (i=0;i<cuts[h].ne;i++){
      for (w=0;w<4;w++)   delete [] edgenod[h][i][w];
      delete [] edgenod[h][i];
    }
    delete [] edgenod[h];
  }
  delete [] edgenod;
  delete [] pom;
  delete [] arsec;
  delete [] cuts;
  delete [] topcuts;

  return 0;
}

