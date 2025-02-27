#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* 
   program vytvori ze souboru se siti ve formatu SIFEL vystupni soubor,
   ktery je potreba pro program partdmesh,ktery rozrezava site na podoblasti
  
   argv[1] - vstupni soubor
   argv[2] - vystupni soubor
   
   190407 JB
*/

int main (int argc, char *argv[])
{
  long i,j,nn,ne,**node,meshtype,a;
  FILE *in,*out;
  char  line [BUFSIZ];
  char *p;
  const char *mezera=" ";
  
  if(argc != 3){
    printf("WARNING: Nesouhlasi pocet paramatru\n");
    printf("vstupni_soubor_se_siti_formatu_jktk   vystupni_soubor\n");
    exit(EXIT_FAILURE);
  }
  printf("Zahajen preprocesing pro partdmesh\n");
  
  in = fopen (argv[1],"r");
  out = fopen (argv[2],"w");
  
  // ************
  // Nacteni uzlu
  // ************
  fscanf (in,"%ld",&nn);
  for(i = 0; i < nn+1; i++){
    fgets(line,BUFSIZ,in); 
  }
  
  // **********************************
  // Nacteni elementu podle jejich typu
  // **********************************
  fscanf (in,"%ld",&ne);
  fprintf (out,"%ld 1\n",ne);
  for (i=0;i<ne;i++){
    fscanf (in,"%ld %ld",&j,&meshtype);
    node = new long* [3];
    fscanf (in,"%ld %ld %ld",node+0,node+1,node+2);
    fscanf (in,"%ld %ld %ld %ld %ld",&j,&j,&j,&j,&j);
    fprintf (out,"%ld %ld %ld\n",node[0],node[1],node[2]);
  }
    
    /*
  switch(meshtype){
    // *********
    // Triangels
    // *********
  case 3:
     printf ("Triangle elements\n");
     fprintf (out,"%ld 1\n",ne);
     printf ("Pocet elementu: %ld\n",ne);
     node = new long* [ne];
     fgets(line,BUFSIZ,in);
     for (i=0;i<ne;i++){
       j = 0;
       node[i] = new long [9];
       fgets(line,BUFSIZ,in);
       // 3 
       p = strtok(line,mezera);
       a = atoi(p);
       node[i][0] = a;
       while (p != NULL){
	 a = atoi(p);
	 node[i][j] = a;
	 j++;
	 p = strtok(NULL,mezera);  
       }  
     }
     for (i=0;i<ne;i++){
       for (j=1;j<4;j++){
	 fprintf (out,"%ld ",node[i][j]);
       }
       fprintf (out,"\n"); 
     }
     break;
     
     // ************** 
     // Quadrilaterals
     // **************
  case 5:
    printf ("Quadrilateral elements\n");
    fprintf (out,"%ld 4\n",ne);
    printf ("Pocet elementu: %ld\n",ne);
    node = new long* [ne];
    fgets(line,BUFSIZ,in);
    for (i=0;i<ne;i++){
      j = 0;
      node[i] = new long [11];
      fgets(line,BUFSIZ,in);
      // 3 
      p = strtok(line,mezera);
      a = atoi(p);
      node[i][0] = a;
      while (p != NULL){
	a = atoi(p);
	node[i][j] = a;
	j++;
	p = strtok(NULL,mezera);  
      }  
    }
    for (i=0;i<ne;i++){
      for (j=1;j<5;j++){
	fprintf (out,"%ld ",node[i][j]);
      }
      fprintf (out,"\n"); 
    }
    break;
    
    // ************
    // Tetrahedrons
    // ************
  case 7:
    printf ("Tetrahedron elements\n");
    fprintf (out,"%ld 2\n",ne);
    printf ("Pocet elementu: %ld\n",ne);
    node = new long* [ne];
    for (i=0;i<ne;i++){
      node[i] = new long [4];
      fscanf (in,"%ld",&j);
      for (j=0;j<4;j++){
	fscanf (in,"%ld",&node[i][j]);
      }
      fscanf (in,"%ld",&j);
    }
    for (i=0;i<ne;i++){
      for (j=0;j<4;j++){
	fprintf (out,"%ld ",node[i][j]);
      }
      fprintf (out,"\n");
    } 
    break;
    
    // ***********
    // Hexahedrons
    // ***********
  case 13:
    printf ("Hexahedron elements\n");
    fprintf (out,"%ld 3\n",ne);
    printf ("Pocet elementu: %ld\n",ne);
    node = new long* [ne];
    for (i=0;i<ne;i++){
      node[i] = new long [8];
      fscanf (in,"%ld",&j);
      for (j=0;j<8;j++){
	fscanf (in,"%ld",&node[i][j]);
      }
      fscanf (in,"%ld",&j);
    }
    for (i=0;i<ne;i++){
      for (j=0;j<8;j++){
	fprintf (out,"%ld ",node[i][j]);
      }
      fprintf (out,"\n");
    }
    break;
  }
    */
  
  printf("Ukoncen preprocesing pro partdmesh\nVytvoren soubor %s se vstupem pro partdmesh\n",argv[2]);
  /*
  // mazani pouzitych poli
  for (i=0;i<ne;i++){
    delete []node[i];
  }
  delete []node;
  */
  fclose (in);
  fclose (out);
}
