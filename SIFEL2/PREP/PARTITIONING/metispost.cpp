# include <stdio.h>
# include <stdlib.h>
# include <string.h>

/* 
   program zpracuje vystup z programu partdmesh 
   vytvori ze souboru se siti ve formatu SIFEL a vystupniho souboru .epart. z partdmesh
   jednotlive soubory *.top s rozrezanou siti ve formatu SIFEL
   
   argv[1] - vstupni soubor s celou nerozrezanou siti
   argv[2] - vystupni soubor z partdmesh *.epart.*
   argv[3] - nd - pocet domen
   argv[4] - spolecne jmeno vystupnich souboru s rozrezanou siti ve formatu SIFEL
   argv[5] - md - mesh description - 1 - all nodes - vsechny uzly maji globalni cislo, to odpovida pubodnimu cislu uzlu
                                         na nerozrezane siti, hranicni uzly maji toto cislo zaporne
                                     2 - boundary nodes - pouze hranicni uzly maji globalni cislo, ostatni maji 0
				    
   190407 JB
*/


int main (int argc, char *argv[])
{
  if( argc != 6){
    printf("Pozor: nesouhlasi pocet vstupnich parametru!\n");
    //printf("soubor se siti pro sekvencni vypocet,spolecny nazev vystupu z metisu,pocet domen,nazev vystupnich souboru\n");
    printf ("soubor s celou siti, vystupni soubor z METISu, pocet podoblasti, nazvy vystupnich souboru bez koncovek, typ site (md)");
    exit(EXIT_FAILURE);
  }
  
  
  long i,j,k,l,n,a,ii,jj;
  long nn,ne,nd,ned,nnd,nbn,meshtype,nel,md;
  long *nprop,*part_elements,*lned,*nod,*multip;
  long **elements;
  double *x,*y,*z;
  FILE *in,*out;
  char *p;
  const char *mezera=" ";
  char  line [BUFSIZ],outt[BUFSIZ],pout[BUFSIZ];
  char pripona2[]=".top";
  
  //  nd - pocet podoblasti
  nd=atoi(argv[3]);
  
  // md - mesh description
  md = atoi(argv[5]);
  
  printf("Zahajen postprocesing vystupu z partdmesh pro %ld domen(y) a ",nd);
  if(md == 1){
    printf("mesh description - all nodes\n");
  }
  else{
    printf("mesh description - boundary nodes\n");
  }
  
  // definice vystupniho souboru - jeho cast
  strcpy(pout,argv[4]);
  
  //
  //  otevreni souboru s celou siti
  //
  in=fopen(argv[1],"rt");
  
  //  nn pocet uzlu cele site
  fscanf(in,"%ld",&nn);
  
  //  alokace poli souradnic uzlu
  x = new double [nn];
  y = new double [nn];
  z = new double [nn];
  
  //  alokace pole nodes_property - obsahuje vlastnost node
  nprop = new long[nn];
  
  //  cteni poli x,y,z a nprop
  for( i = 0; i<nn ; i++){
    fscanf(in,"%lf %lf %lf %ld",x+i,y+i,z+i,nprop+i);
  }
  printf("Uzly nacteny ze souboru %s\n",argv[1]);
  
  
  fscanf (in,"%ld %ld",&ne,&meshtype);
  switch(meshtype){
    // triangels
  case 3:
    printf ("Triangel elements\n");
    elements = new long* [ne];
    fgets(line,BUFSIZ,in);
    for (i=0;i<ne;i++){
      j = 0;
      elements[i] = new long [8];
      fgets(line,BUFSIZ,in);
      // 3 
      p = strtok(line,mezera);
      a = atoi(p);
      elements[i][j] = a;
      j++;
      //printf ("%ld  ",elements[i][j]);
      p = strtok(NULL,mezera);  
      while (p != NULL){
	a = atoi(p);
	elements[i][j] = a;
	//printf ("%ld  ",elements[i][j]);
	j++;
	p = strtok(NULL,mezera);  
      }
     //printf ("\n");
    }
    nel = j;
    //printf ("%ld \n",nel);
    break;
    
    // quadrilaterals
  case 5:
    printf ("Quadrilateral elements\n");
    elements = new long* [ne];
    fgets(line,BUFSIZ,in);
    for (i=0;i<ne;i++){
      j = 0;
      elements[i] = new long [11];
      fgets(line,BUFSIZ,in);
      // 3 
      p = strtok(line,mezera);
      a = atoi(p);
      elements[i][j] = a;
      j++;
      p = strtok(NULL,mezera);  
      elements[i][0] = a;
      while (p != NULL){
	a = atoi(p);
	elements[i][j] = a;
	j++;
	p = strtok(NULL,mezera);  
      }  
    }
    nel = j;
    break;
    
    // tetrahedrons
  case 7:
    printf ("Tetrahedron elements\n");
    elements = new long* [ne];
    for (i=0;i<ne;i++){
      elements[i] = new long [6];
      for (j=0;j<6;j++){
	fscanf (in,"%ld",&elements[i][j]);
	//printf("%ld   ",elements[i][j]);
      }
      //printf("\n");
    }
    
    break;
    
    // hexahedrons
  case 13:
    printf ("Hexahedron elements\n");
    elements = new long* [ne];
    for (i=0;i<ne;i++){
      elements[i] = new long [10];
      for (j=0;j<10;j++){
	fscanf (in,"%ld",&elements[i][j]);  
      }
    }
    break;
  }
  
  fclose(in);
  printf("Elementy nacteny ze souboru %s\n",argv[1]);
  
  
  part_elements = new long[ne];
  
  //  cteni vystupniho souboru z METISu
  in=fopen(argv[2],"rt");
  for( i =  0; i<ne ; i++)
    fscanf(in,"%ld",&part_elements[i]);
  fclose(in); 
  printf("Deleni elementu na domeny nacteno ze souboru %s\n",argv[2]);
  
  
  
  
  
  //  ned - pocet elementu na jednotlivych domenach
  printf("Pocet elementu na jednotlivych domenach\n"); 
  lned = new long [nd];
  
  for(i = 0; i < nd; i++){
    lned[i] = 0;
  }
  for(j = 0; j < ne; j++){
    lned[part_elements[j]]++;
  }
  for (i=0;i<nd;i++){
    printf("Domena %ld : %ld\n",i,lned[i]);
  }
  
  //kontrola celkoveho poctu elementu
  n = 0;
  for(i = 0; i < nd; i++) n+=lned[i];
  if (n != ne)
    printf("Soucet prvku na podoblastech neni roven poctu prvku v puvodni siti!\n");
  else
    printf("Soucet prvku na podoblastech je roven poctu prvku v puvodni siti!\n"); 
  
  //  uzlova multiplicita - pocet podoblasti, na kterych lezi dany uzel
  multip = new long [nn];
  for(i = 0; i < nn; i++)
    multip[i]=0;
  
  nod = new long [nn];
  
  switch(meshtype){
    // triangels
  case 3:
    for (i=0;i<nd;i++){
      
      for (j=0;j<nn;j++){
	nod[j]=0;
      }
      
      //  pocet prvku na i-te podoblasti
      ned=0;
      for (j=0;j<ne;j++){
	if (part_elements[j]==i){
	  ned++;
	  ii=part_elements[j];
	  for (k=1;k<4;k++){
	    jj=elements[j][k]-1;
	    nod[jj]++;
	  }
	}
      }
      
      //  pocet uzlu na i-te podoblasti
      nnd=0;
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  nnd++;
	  multip[j]++;
	}
      }
    }
    // mesh description
    switch(md){
      // all nodes
    case 1:
       nbn=1;
       for (i=0;i<nn;i++){
	 if (multip[i]>1){
	   multip[i]=(i+1)*(-1); 
	   nbn++;
	 }
	 else{
	   multip[i]=i+1; 
	 }
       }
      break;
      // boundary nodes
    case 2:
      nbn=1;
      for (i=0;i<nn;i++){
	if (multip[i]>1){
	  multip[i]=nbn;
	  nbn++;
	}
	else{
	  multip[i]=0;
	}
      }
      break;
    }
    
    printf ("Pocet hranicnich uzlu:   %ld\n",nbn);
    
    
    
    // pocet uzlu na jednotlivych domenach
    
    
    for (i=0;i<nd;i++){
      //sprintf(outt,"%s%d",pout,i+1);
      sprintf(outt,"%s.top.%ld",pout,i);
      //strcat(outt,pripona2);	
      printf("Zapisuji vystup do souboru %s\n",outt);
      out=fopen(outt,"wt");
      
      
      for (j=0;j<nn;j++){
	nod[j]=0;
      }
      
      //  pocet prvku na i-te podoblasti
      ned=0;
      for (j=0;j<ne;j++){
	if (part_elements[j]==i){
	  ned++;
	  ii=part_elements[j];
	  for (k=1;k<4;k++){
	    jj=elements[j][k]-1;
	    nod[jj]++;
	  }
	}
      }
      
      //  pocet uzlu na i-te podoblasti
      nnd=0;
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  nnd++;
	}
      }
      
      
      fprintf (out,"%ld\n",nnd);
      k=1;
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  fprintf (out,"%5le   %5le   %5le      %2ld\n",x[j],y[j],z[j],nprop[j]);
	  nod[j]=k;
	  k++;
	}
      }
      
      fprintf (out,"%ld %ld\n",ned,meshtype);
      for (j=0;j<ne;j++){
	if (part_elements[j]==i){
	  fprintf (out,"%5ld",elements[j][0]);
	  for (k=1;k<4;k++){
	    fprintf (out,"   %5ld",nod[elements[j][k]-1]);
	  }
	  if(nel == 9){
	    fprintf (out,"     %5ld   %5ld   %5ld    %5ld    %5ld\n",elements[j][4],elements[j][5],elements[j][6],elements[j][7],elements[j][8]);
	  }
	  else{
	    fprintf (out,"     %5ld\n",elements[j][4]);
	  }
	}
      }
      
      fprintf (out,"\n");
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  fprintf (out,"%5ld %5ld\n",multip[j]);
	}
      }
      
      
      fclose (out);
    }
    break;
    // quadrilateral
  case 5:
    for (i=0;i<nd;i++){
      
      for (j=0;j<nn;j++){
	nod[j]=0;
      }
      
      //  pocet prvku na i-te podoblasti
      ned=0;
      for (j=0;j<ne;j++){
	if (part_elements[j]==i){
	  ned++;
	  ii=part_elements[j];
	  for (k=1;k<5;k++){
	    jj=elements[j][k]-1;
	    nod[jj]++;
	  }
	}
      }
      
      //  pocet uzlu na i-te podoblasti
      nnd=0;
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  nnd++;
	  multip[j]++;
	}
      }
    }
    // mesh description
    switch(md){
      // all nodes
    case 1:
       nbn=1;
       for (i=0;i<nn;i++){
	 if (multip[i]>1){
	   multip[i]=(i+1)*(-1); 
	   nbn++;
	 }
	 else{
	   multip[i]=i+1; 
	 }
       }
      break;
      // boundary nodes
    case 2:
      nbn=1;
      for (i=0;i<nn;i++){
	if (multip[i]>1){
	  multip[i]=nbn;
	  nbn++;
	}
	else{
	  multip[i]=0;
	}
      }
      break;
    }
    
    printf ("Pocet hranicnich uzlu:   %ld\n",nbn);
    
    
    
    // pocet uzlu na jednotlivych domenach
    
    
    for (i=0;i<nd;i++){
      //sprintf(outt,"%s%d",pout,i+1);
      //strcat(outt,pripona2);
      sprintf(outt,"%s.top.%ld",pout,i);
      
      printf("Zapisuji vystup do souboru %s\n",outt);
      out=fopen(outt,"wt");
      
      
      for (j=0;j<nn;j++){
	nod[j]=0;
      }
      
      //  pocet prvku na i-te podoblasti
      ned=0;
      for (j=0;j<ne;j++){
	if (part_elements[j]==i){
	  ned++;
	  ii=part_elements[j];
	  for (k=1;k<5;k++){
	    jj=elements[j][k]-1;
	    nod[jj]++;
	  }
	}
      }
      
      //  pocet uzlu na i-te podoblasti
      nnd=0;
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  nnd++;
	}
      }
      
      
      fprintf (out,"%ld\n",nnd);
      k=1;
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  fprintf (out,"%5le   %5le   %5le      %2ld\n",x[j],y[j],z[j],nprop[j]);
	  nod[j]=k;
	  k++;
	}
      }
      
      fprintf (out,"%ld %ld\n",ned,meshtype);
      for (j=0;j<ne;j++){
	if (part_elements[j]==i){
	  fprintf (out,"%5ld",elements[j][0]);
	  for (k=1;k<5;k++){
	    fprintf (out,"   %5ld",nod[elements[j][k]-1]);
	  }
	  if(nel == 11){
	    fprintf (out,"     %5ld    %5ld    %5ld    %5ld   %5ld\n",elements[j][5],elements[j][6],elements[j][7],elements[j][8],elements[j][9] ,elements[j][10]);
	  }
	  else{
	    fprintf (out,"     %5ld\n",elements[j][5]);
	  }
	}
      }
      
      fprintf (out,"\n");
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  fprintf (out,"%5ld\n",multip[j]);
	}
      }
      
      
      fclose (out);
    }
    
    break;
    
    // tetrahedron
  case 7:
    for (i=0;i<nd;i++){
      
      for (j=0;j<nn;j++){
	nod[j]=0;
      }
      
      //  pocet prvku na i-te podoblasti
      ned=0;
      for (j=0;j<ne;j++){
	if (part_elements[j]==i){
	  ned++;
	  ii=part_elements[j];
	  for (k=1;k<5;k++){
	    jj=elements[j][k]-1;
	    nod[jj]++;
	  }
	}
      }
      
      //  pocet uzlu na i-te podoblasti
      nnd=0;
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  nnd++;
	  multip[j]++;
	}
      }
    }
    
     // mesh description
    switch(md){
      // all nodes
    case 1:
       nbn=1;
       for (i=0;i<nn;i++){
	 if (multip[i]>1){
	   multip[i]=(i+1)*(-1); 
	   nbn++;
	 }
	 else{
	   multip[i]=i+1; 
	 }
       }
      break;
      // boundary nodes
    case 2:
      nbn=1;
      for (i=0;i<nn;i++){
	if (multip[i]>1){
	  multip[i]=nbn;
	  nbn++;
	}
	else{
	  multip[i]=0;
	}
      }
      break;
    }
        
    printf ("Pocet hranicnich uzlu   %ld\n",nbn);
    
    
    
    // pocet uzlu na jednotlivych domenach
    
    //long *lnnd;
    
    for (i=0;i<nd;i++){
      //sprintf(outt,"%s%d",pout,i+1);
      //strcat(outt,pripona2);	
      sprintf(outt,"%s.top.%ld",pout,i);
      printf("Zapisuji vystup do souboru %s\n",outt);
      out=fopen(outt,"wt");
      
      
      for (j=0;j<nn;j++){
	nod[j]=0;
      }
      
      //  pocet prvku na i-te podoblasti
      ned=0;
      for (j=0;j<ne;j++){
	if (part_elements[j]==i){
	  ned++;
	  ii=part_elements[j];
	  for (k=1;k<5;k++){
	    jj=elements[j][k]-1;
	    nod[jj]++;
	  }
	}
      }
      
      //  pocet uzlu na i-te podoblasti
      nnd=0;
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  nnd++;
	}
      }
      
      
      fprintf (out,"%ld\n",nnd);
      k=1;
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  fprintf (out,"%5le   %5le   %5le      %2ld\n",x[j],y[j],z[j],nprop[j]);
	  nod[j]=k;
	  k++;
	}
      }
      
      fprintf (out,"%ld %ld\n",ned,meshtype);
      for (j=0;j<ne;j++){
	if (part_elements[j]==i){
	  fprintf (out,"%5ld",elements[j][0]);
	  for (k=1;k<5;k++){
	    fprintf (out,"   %5ld",nod[elements[j][k]-1]);
	  }
	  fprintf (out,"     %5ld\n",elements[j][5]);
	}
      }
      
      fprintf (out,"\n");
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  fprintf (out,"%5ld\n",multip[j]);
	}
      }
      
      
      fclose (out);
    }
    break;
    
    // hexahedron
  case 13:
    for (i=0;i<nd;i++){
      
      for (j=0;j<nn;j++){
	nod[j]=0;
      }
      
      //  pocet prvku na i-te podoblasti
      ned=0;
      for (j=0;j<ne;j++){
	if (part_elements[j]==i){
	  ned++;
	  ii=part_elements[j];
	  for (k=1;k<9;k++){
	    jj=elements[j][k]-1;
	    nod[jj]++;
	  }
	}
      }
      
      //  pocet uzlu na i-te podoblasti
      nnd=0;
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  nnd++;
	  multip[j]++;
	}
      }
    }
    
     // mesh description
    switch(md){
      // all nodes
    case 1:
       nbn=1;
       for (i=0;i<nn;i++){
	 if (multip[i]>1){
	   multip[i]=(i+1)*(-1); 
	   nbn++;
	 }
	 else{
	   multip[i]=i+1; 
	 }
       }
      break;
      // boundary nodes
    case 2:
      nbn=1;
      for (i=0;i<nn;i++){
	if (multip[i]>1){
	  multip[i]=nbn;
	  nbn++;
	}
	else{
	  multip[i]=0;
	}
      }
      break;
    }
        
    printf ("Pocet hranicnich uzlu   %ld\n",nbn);
    
    
    
    // pocet uzlu na jednotlivych domenach
    
    //long *lnnd;
    
    for (i=0;i<nd;i++){
      //sprintf(outt,"%s%d",pout,i+1);
      //strcat(outt,pripona2);	
      sprintf(outt,"%s.top.%ld",pout,i);
      printf("Zapisuji vystup do souboru %s\n",outt);
      out=fopen(outt,"wt");
      
      
      for (j=0;j<nn;j++){
	nod[j]=0;
      }
      
      //  pocet prvku na i-te podoblasti
      ned=0;
      for (j=0;j<ne;j++){
	if (part_elements[j]==i){
	  ned++;
	  ii=part_elements[j];
	  for (k=1;k<9;k++){
	    jj=elements[j][k]-1;
	    nod[jj]++;
	  }
	}
      }
      
      //  pocet uzlu na i-te podoblasti
      nnd=0;
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  nnd++;
	}
      }
      
      
      fprintf (out,"%ld\n",nnd);
      k=1;
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  fprintf (out,"%5le   %5le   %5le      %2ld\n",x[j],y[j],z[j],nprop[j]);
	  nod[j]=k;
	  k++;
	}
      }
      
      fprintf (out,"%ld %ld\n",ned,meshtype);
      for (j=0;j<ne;j++){
	if (part_elements[j]==i){
	  fprintf (out,"%5ld",elements[j][0]);
	  for (k=1;k<9;k++){
	    fprintf (out,"   %5ld",nod[elements[j][k]-1]);
	  }
	  fprintf (out,"     %5ld\n",elements[j][9]);
	}
      }
      
      fprintf (out,"\n");
      for (j=0;j<nn;j++){
	if (nod[j]>0){
	  fprintf (out,"%5ld\n",multip[j]);
	}
      }
      
      
      fclose (out);
    }
    break;
  }
  printf("Postprocesing vystupu z partdmesh ukoncen\n");
  
  // mazani pouzitych poli
  delete []nprop;
  delete []part_elements;
  delete []lned;
  delete []nod;
  delete []multip;
  for(i = 0; i < ne; i++) delete []elements[i];
  delete []elements;
  delete []x;
  delete []y;
  delete []z;
  


}




  
