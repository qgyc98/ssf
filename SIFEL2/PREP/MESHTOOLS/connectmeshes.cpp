#include "connectmeshes.h"

void check_reqsec(XFILE *in)
 {
   long err;

   // index of sections must be created
   if (in->index_created != 1)
     {
       print_err("Index of sections has not beeen created", __FILE__, __LINE__, __func__);
       abort();
     }
   
   // section with input files must be detected
   err = xf_setsec(in, bsec_str[begsec_files]);  
   switch (err)
       {
       case 1:
      print_err("Section with input files has not been detected", __FILE__, __LINE__, __func__);
      abort();
     default:
      break;
  }
}


connectmeshes::connectmeshes(int argc,char *argv[])
{
  Argc = argc;
  Argv = argv;  

}
connectmeshes::~connectmeshes()
{
  long i;
  // for(i = 0; i < Argc; i++){
//     delete []Argv[i];
//   }
//   delete []Argv;
  
  for(i = 0; i < nfiles; i++){
    delete []interfnodes[i];
  }
  delete []ltg;
  delete []interfnodes;
  
}

int connectmeshes::read_input_data()
{
  printf("Reading of input data parametrs\n");
  long l;
  in = xfopen(Argv[1],"r");
  if (in == NULL){
    fprintf (stderr,"\n Input file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    delete Gtm; delete Top;
    return(2);
  }
  logname = new char[strlen(Argv[1])+5];
  strcpy(logname, Argv[1]);
  l = strlen(logname);
  tmp = strrchr(logname, '.');
  if (tmp)
    sprintf(tmp+1, "plg");    
  else
    sprintf(logname+l, "plg");
  Log = fopen(logname, "w");
  delete [] logname;
  in->warning = 1;

  in->kwdmode = sect_mode_seq;
  //in->kwdmode = sequent_mode;
  in->ignorecase = 1;
  // detection of sections in the preprocesor file
  xfdetect_sect(in, bsec_kwdset, esec_kwdset);
  // checking of sections that have to be detected
  check_reqsec(in);

  // reading of input topology 
  xf_setsec(in, bsec_str[begsec_files]);
  // reading of line with topology file name
  xfscanf(in, "%k%s","topology_file",&topfile);
  // reading of line with topology file format indicator,
  xfscanf(in, "%k%m", "mesh_format", &meshform_kwdset, &d.meshfmt);
  // reading of line with edge number indicator
  xfscanf(in, "%k%ld", "edge_numbering", &d.redgn);
  // reading of line with edge number indicator
  // reading of line with output file name
  xfscanf(in,"%k%s","output_file_name",&outputfname);
  // reading of tolerance for coordinate difference under which the nodes are considered to be identical
  xfscanf(in, "%k%le", "tolerance", &coordtol);
  
  //  this is the sequential version
  d.paral=0;
  // reading of line with edge number indicator
  xfscanf(in, "%k%ld", "number_of_files", &nfiles);
  
  xfclose(in);
  //printf(" . . . OK\n");
  return(0);
}

  
/**
  Function reads topology file in, format of the file is described by structure d.

  @param in - pointer to the opened XFILE structure
  @param d  - structure with input data format description


  Function returns:
  @retval 0 : on success
  @retval 4 : in case of unknown mesh format

  In case of reading alternative file format :
  @retval 1 : on error in reading of node
  @retval 2 : on error in reading of element
  @retval 3 : on error in reading of global node numbers
*/

int connectmeshes::read_topology(long nfil)
{
  double ts,te;
  
  // set tprob for reading of mesh
  ts = clock();
  topf = xfopen(d.topf,"r");
  if (topf == NULL){
    fprintf (stderr,"\n Topology file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(3);
  }
  topf->warning = 1;
  topf->kwdmode = ignore_kwd;
  topf->ignorecase = 1;
  err = input_siftop(nfil);
  xfclose(topf);
  if (err){
    print_err("\nReading of mesh topology failed\n", __FILE__, __LINE__, __func__);
    return(3);
  }
  else{
    //printf(" OK\n");
  }
  te = clock();
  //printf("Time of reading topology %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  return(0);
  
  
}

long connectmeshes::input_siftop(long nfil)
{
  long ret;
  
  //fprintf(stdout, "Reading of mesh topology . . .");
  switch (d.meshfmt)
  {
    case t3d:
      SifTop[nfil].import_t3d(topf, d.paral);
      break;
    case sifel:
      ret = SifTop[nfil].read(topf, d.paral, d.redgn);
      return(ret);
    default:
      print_err("unknown mesh format is required", __FILE__, __LINE__, __func__);
      return(4);
  }
  return(0);
}

void connectmeshes::establish_ltg()
{
  //printf("--------------------------\n");
  fprintf(stdout, "Searching interface nodes\n");
  long i,j,k,l,m,n,p;
  long nedges,nsurf;
  long *edgenn,*dupledges,*aux,*duplsurf,*surfnn;
  long **edgenodes,**surfnodes;
  long auxninterfnodes;
  long *auxinterfnodes;
  long *duplnodes;
  long *elemedge;
  long *elemsurf;
  double ts,te,tss,tee;
  long pointer1,pointer2;
  long ns1,ns2;
    
  ts = clock();
  n = 0;
  for(i = 0; i < nfiles; i++){
    n+=Gtop[i].nn;
  }
  ltg=new long[n];
  n = 0;
  for(i = 0; i < nfiles; i++){
    for(j = 0; j < Gtop[i].nn; j++){
      ltg[n] = 0;
      n++;
    }
  }
  
  ninterfnodes = new long[nfiles];
  interfnodes = new long*[nfiles];
  for(i = 0; i < nfiles; i++){
    // number of all edges in mesh
    nedges = 0;
    nsurf = 0;
    elemedge = new long[Gtop[i].ne];
    elemsurf = new long[Gtop[i].ne];
    for(j = 0; j < Gtop[i].ne; j++){
      switch(Gtop[i].gelements[j].get){
      case linbar:
      case quadbar:{
	break;
      }
      case lintriag:
      case quadtriag:
      case linquad:
      case quadquad:{
	// number of edge on i-th element - gtopology
	elemedge[j] = nedges;
	nedges+= Gtop[i].give_ned(j);
	break;
      }
      case lintetra:
      case quadtetra:
      case linhexa:
      case quadhexa:{
	// number of edge on i-th element - gtopology
	elemsurf[j] = nsurf;
	nsurf+= Gtop[i].give_nsurf(j);
	break;
      }
      case noelem:
      default:{
	print_err("Unknown element type", __FILE__, __LINE__, __func__);
	break;
      }
      }
    }
    //fprintf(Log,"pocet hran %ld pocet ploch %ld\n",nedges,nsurf);
    
    //     tss = clock();
    // pomocne pole
    aux = new long[8];
    if(nedges != 0){
      edgenn = new long[nedges];
      edgenodes = new long*[nedges];
      dupledges = new long[nedges];
    }
    if(nsurf != 0){
      surfnn = new long[nsurf];
      surfnodes = new long*[nsurf];
      duplsurf = new long[nsurf];
    }
    l = 0;
    n = 0;
    for(j = 0; j < Gtop[i].ne; j++){
      switch(Gtop[i].gelements[j].get){
      case  linbar:
      case quadbar:{
	break;
      }
      case lintriag:
      case quadtriag:
      case linquad:
      case quadquad:{
	for(k = 0; k < Gtop[i].give_ned(j); k++ ){
	  edgenn[l] = Gtop[i].give_nned(j);
	  edgenodes[l] = new long[edgenn[l]];
	  Gtop[i].give_edge_nodes(j,k,edgenodes[l]);
	  // trideni pole s uzly na hrane(od nejmensiho po nejvetsi)
	  if(edgenn[l] == 2){
	    if(edgenodes[l][0] > edgenodes[l][1]){
	      m = edgenodes[l][1];
	      edgenodes[l][1] = edgenodes[l][0];
	      edgenodes[l][0] = m;
	    }
	  }
	  if(edgenn[l] == 3){
	    if(edgenodes[l][0] > edgenodes[l][2]){
	      m = edgenodes[l][2];
	      edgenodes[l][2] = edgenodes[l][0];
	      edgenodes[l][0] = m;
	    }
	    if(edgenodes[l][0] > edgenodes[l][1]){
	      m = edgenodes[k][1];
	      edgenodes[l][1] = edgenodes[l][0];
	      edgenodes[l][0] = m;
	    }
	    if(edgenodes[l][1] > edgenodes[l][2]){
	      m = edgenodes[l][2];
	      edgenodes[l][2] = edgenodes[l][1];
	      edgenodes[l][1] = m;
	    }
	  }
	  dupledges[l] = 0;
	  l++;
	}
	break;
      }
      case lintetra:
      case quadtetra:
      case linhexa:
      case quadhexa:{
	for(k = 0; k < Gtop[i].give_nsurf(j); k++ ){
	  surfnn[n] = Gtop[i].give_nnsurf(j);
	  surfnodes[n] = new long[surfnn[n]];
	  Gtop[i].give_surf_nodes(j,k,surfnodes[n]);
	  // bubble sorting
	  for(l =  Gtop[i].give_nnsurf(j) - 1; l > 0; l--){
	    for(p = 0; p < l; p++){
	      if( surfnodes[n][p] >  surfnodes[n][p+1]){
		m = surfnodes[n][p];   
		surfnodes[n][p] = surfnodes[n][p+1];
		surfnodes[n][p+1] = m;
	      }
	    }
	  }
	  duplsurf[n] = 0;
	  n++;
	}
	break;
      }
      case noelem:
      default:{
	print_err("Unknown element type", __FILE__, __LINE__, __func__);
	break;
      }
      }
    }
    //tee = clock();
    //printf("List of nodes on edges or/and surfaces is created\n",(tee-tss)/(double)CLOCKS_PER_SEC);
    //printf("List of nodes on edges or/and surfaces is created\n");
    tss = clock();
    for(j = 0; j < Gtop[i].ne; j++){
      switch(Gtop[i].gelements[j].get){
      case linbar:
      case quadbar:{
	break;
      }
      case lintriag:
      case quadtriag:
      case linquad:
      case quadquad:{
	ns1=Gtop[i].give_ned(j);
	pointer1=elemedge[j];
	for(k = 0; k < Gtop[i].nadjelel[j]; k++){
	  if(Gtop[i].adjelel[j][k] != j ){
	    pointer2=elemedge[Gtop[i].adjelel[j][k]];
	    ns2=Gtop[i].give_ned(Gtop[i].adjelel[j][k]);
	    for(l = 0; l < ns1; l++){
	      for(m = 0; m < ns2; m++){
		if(edgenodes[pointer1+l][0] == edgenodes[pointer2+m][0]){
		  if(edgenn[pointer1+l] == 2 && edgenn[pointer2+m] == 2){
		    if(edgenodes[pointer1+l][1] == edgenodes[pointer2+m][1]){
		      dupledges[pointer1+l]++;
		      dupledges[pointer2+m]++;
		      break;
		    }
		  }
		  if(edgenn[pointer1+l] == 3 && edgenn[pointer2+m] == 3){
		    if(edgenodes[pointer1+l][1] == edgenodes[pointer2+m][1]){
		      if(edgenodes[pointer1+l][2] == edgenodes[pointer2+m][2]){
			dupledges[pointer1+l]++;
			dupledges[pointer2+m]++;
			break;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	break;
      }
      case lintetra:
      case quadtetra:
      case linhexa:
      case quadhexa:{
	ns1=Gtop[i].give_nsurf(j);
 	pointer1=elemsurf[j];
	for(k = 0; k < Gtop[i].nadjelel[j];k++){
 	  if(Gtop[i].adjelel[j][k] != j ){
	    pointer2=elemsurf[Gtop[i].adjelel[j][k]];
 	    ns2=Gtop[i].give_nsurf(Gtop[i].adjelel[j][k]);
	    for(l = 0; l < ns1; l++){
 	      for(m = 0; m < ns2; m++){
		if(surfnodes[pointer1+l][0] == surfnodes[pointer2+m][0]){
		  if(surfnn[pointer1+l] == 3 && surfnn[pointer2+m] == 3){
		    if(surfnodes[pointer1+l][1] == surfnodes[pointer2+m][1]){
		      if(surfnodes[pointer1+l][2] == surfnodes[pointer2+m][2]){
			duplsurf[pointer1+l]++;
			duplsurf[pointer2+m]++;
			break;
		      }
		    }
		  }
		  if(surfnn[pointer1+l] == 4 && surfnn[pointer2+m] == 4){
		    if(surfnodes[pointer1+l][1] == surfnodes[pointer2+m][1]){
		      if(surfnodes[pointer1+l][2] == surfnodes[pointer2+m][2]){
			if(surfnodes[pointer1+l][3] == surfnodes[pointer2+m][3]){
			  duplsurf[pointer1+l]++;
			  duplsurf[pointer2+m]++;
			  break;
			}
		      }
		    }
		  }
		  if(surfnn[pointer1+l] == 6 && surfnn[pointer2+m] == 6){
		    if(surfnodes[pointer1+l][1] == surfnodes[pointer2+m][1]){
		      if(surfnodes[pointer1+l][2] == surfnodes[pointer2+m][2]){
			if(surfnodes[pointer1+l][3] == surfnodes[pointer2+m][3]){
			  if(surfnodes[pointer1+l][4] == surfnodes[pointer2+m][4]){
			    if(surfnodes[pointer1+l][5] == surfnodes[pointer2+m][5]){
			      duplsurf[pointer1+l]++;
			      duplsurf[pointer2+m]++;
			      break;
			    }
			  }
			}
		      }
		    }
		  }
		  if(surfnn[pointer1+l] == 8 && surfnn[pointer2+m] == 8){
		    if(surfnodes[pointer1+l][1] == surfnodes[pointer2+m][1]){
		      if(surfnodes[pointer1+l][2] == surfnodes[pointer2+m][2]){
			if(surfnodes[pointer1+l][3] == surfnodes[pointer2+m][3]){
			  if(surfnodes[pointer1+l][4] == surfnodes[pointer2+m][4]){
			    if(surfnodes[pointer1+l][5] == surfnodes[pointer2+m][5]){
			      if(surfnodes[pointer1+l][6] == surfnodes[pointer2+m][6]){
				if(surfnodes[pointer1+l][7] == surfnodes[pointer2+m][7]){
				  duplsurf[pointer1+l]++;
				  duplsurf[pointer2+m]++;
				  break;
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
 	      }
 	    }
 	  }
 	}
	pointer1+=Gtop[i].give_nnsurf(j);
	break;
      }
      case noelem:
      default:{
	print_err("Unknown element type", __FILE__, __LINE__, __func__);
	break;
      }
      }
    }
    
    
    tee = clock();
    //printf("List of potencial contact nodes creates %le s\n",(tee-tss)/(double)CLOCKS_PER_SEC);
    //     printf("List of potencial contact nodes is created\n");
    tss = clock();
    auxninterfnodes = 0;
    if(nedges != 0){
      //   for(j = 0; j < nedges; j++){
      // 	fprintf(Log,"edge %ld %ld\n",j,dupledges[j]);
      //       }
      for(j = 0; j < nedges; j++){
	if(dupledges[j] == 0){
	  auxninterfnodes+=edgenn[j];
	}
      }
    }
    if(nsurf != 0){
      for(j = 0; j < nsurf; j++){
	if(duplsurf[j] == 0){
	  auxninterfnodes+=surfnn[j];
	}
      }
    }
    //fprintf(Log,"auxninterfnodes %ld\n",auxninterfnodes);
    m = 0;
    auxinterfnodes = new long[auxninterfnodes];
    if(nedges != 0){
      for(j = 0; j < nedges; j++){
	if(dupledges[j] == 0){
	  for(k = 0; k < edgenn[j]; k++){
	    auxinterfnodes[m] = edgenodes[j][k];
	    m++;
	  }
	}
      }
      delete []dupledges;
      delete []edgenn;
      for(j = 0; j < nedges; j++){
	delete []edgenodes[j];
      }
      delete []edgenodes;
    }
    if(nsurf != 0){
      for(j = 0; j < nsurf; j++){
	if(duplsurf[j] == 0){
	  for(k = 0; k < surfnn[j]; k++){
	    auxinterfnodes[m] = surfnodes[j][k];
	    m++;
	  }
	}
      }
      delete []duplsurf;
      delete []surfnn;
      for(j = 0; j < nsurf; j++){
	delete []surfnodes[j];
      }
      delete []surfnodes;
      
    }
     for(j = 0; j < auxninterfnodes; j++){
           fprintf(Log,"int nodes %ld\n",auxinterfnodes[j]);
         }
    duplnodes = new long [auxninterfnodes];
    for(j = 0; j < auxninterfnodes; j++){
      duplnodes[j] = 0;
    }
    
    // hledani dulpikovanych uzlu
    for(j = 0; j < auxninterfnodes; j++){
      for(k = j+1; k < auxninterfnodes; k++){
	if(auxinterfnodes[j] == auxinterfnodes[k]){
	  duplnodes[k]++;
	  break;
	}
      }
    }
    tee = clock();
    printf("hledani duplicitnich uzlu %le s\n",(tee-tss)/(double)CLOCKS_PER_SEC);
    
    //kontrolni tisk
    for(j = 0; j < auxninterfnodes; j++){
        fprintf(Log,"%ld duplnodes %ld\n",j,duplnodes[j]);
    }
    
    tss = clock();
    ninterfnodes[i] = 0;
    for(j = 0; j < auxninterfnodes; j++){
      if( duplnodes[j] == 0){
	ninterfnodes[i]++;
      }
    }
    m = 0;
    interfnodes[i] = new long[ninterfnodes[i]];
    for(j = 0; j < auxninterfnodes; j++){
      if(duplnodes[j] == 0){
	interfnodes[i][m]=auxinterfnodes[j];
	m++;
      }
    }
    delete []duplnodes;
    delete []auxinterfnodes;
  }
  tee = clock();
  printf("pocet hranicnich %ld\n",m);
  printf("seznam hranicnich %le s\n",(tee-tss)/(double)CLOCKS_PER_SEC);
  //  printf("seznam hranicnich %le s\n",(tee-tss)/(double)CLOCKS_PER_SEC);
   for(i = 0 ; i < nfiles; i++){
     fprintf(Log,"#  interface nodes %ld\n",ninterfnodes[i]);
     for(j = 0 ; j < ninterfnodes[i]; j++){
       fprintf(Log,"%ld\n",interfnodes[i][j]);
     }
   }

  tss = clock();
  m = -1;
  p = 0;
  long u,v;
  for(i = 0; i < nfiles; i++){
    for(k = 0; k < ninterfnodes[i]; k++){
      u = p+Gtop[i].nn;
      for(j = i+1; j < nfiles; j++){
	for(l = 0; l < ninterfnodes[j]; l++){
	  if(fabs(Gtop[i].gnodes[interfnodes[i][k]].x - Gtop[j].gnodes[interfnodes[j][l]].x) <= coordtol){
	    if(fabs(Gtop[i].gnodes[interfnodes[i][k]].y - Gtop[j].gnodes[interfnodes[j][l]].y) <= coordtol){
	      if(fabs(Gtop[i].gnodes[interfnodes[i][k]].z - Gtop[j].gnodes[interfnodes[j][l]].z) <= coordtol){
		n=p+interfnodes[i][k];
		v=u+interfnodes[j][l];
		if(ltg[n]== 0){
		  ltg[n] = m;
		  m--;
		}
		if(ltg[v]== 0){
   		  ltg[v]= n;
   		}
  		printf("shoda %ld %ld %ld %ld\n",interfnodes[i][k],interfnodes[j][l],n,v);
		break;
	      }
	    }
	  }
	}
	u+=Gtop[j].nn;
      }
      n++;
    }
    p+=Gtop[i].nn;
  }
  //printf("m %ld\n",m);
  // kontrolni tisk
   n = 0;
   for(i = 0; i < nfiles; i++){
     fprintf(Log,"domena %ld\n",i);
     for(j = 0; j < SifTop[i].nn; j++){
       fprintf(Log,"ltg[%ld] = %ld\n",n,ltg[n]);
       n++;
     }
  }
  tee = clock();
  //printf("ltg %le s\n",(tee-tss)/(double)CLOCKS_PER_SEC);
  n = 0;
  for(i = 0; i < nfiles; i++){
    for(j = 0; j < SifTop[i].nn; j++){
      if(ltg[n] > 0){
	ltg[n]=ltg[ltg[n]];
      }
    n++;
    }
  }
  
  n = 0;
  m = 0;
  for(i = 0; i < nfiles; i++){
    for(j = 0; j < SifTop[i].nn; j++){
      if(ltg[n] == 0){
	ltg[n]=m;
	m++;
      }
    n++;
    }
  }

//   n = 0;
//   for(i = 0; i < nfiles; i++){
//     fprintf(Log,"domena %ld\n",i);
//     for(j = 0; j < SifTop[i].nn; j++){
//       fprintf(Log,"ltg[%ld] = %ld\n",n,ltg[n]);
//       n++;
//     }
//   }
  
  n = 0;
  for(i = 0; i < nfiles; i++){
    for(j = 0; j < SifTop[i].nn; j++){
      if(ltg[n] < 0){
	ltg[n]=m+(-1)*ltg[n]-1;
      }
      n++;
    }
  }
  
//   n = 0;
//   for(i = 0; i < nfiles; i++){
//     fprintf(Log,"domena %ld\n",i);
//     for(j = 0; j < SifTop[i].nn; j++){
//       fprintf(Log,"ltg[%ld] = %ld\n",n,ltg[n]);
//       n++;
//     }
//   }
    
  
  
  delete []Gtop;
}


void connectmeshes::print()
{
  //printf("--------------------------\n");
  printf("File %s with output mesh was printed\n",outputfname);
  
  
  
  long i, j;
  long renn,rene,nn,ne;
  FILE *out;
  
  out=fopen(outputfname,"w");

  nn = 0;
  for(i = 0; i < nfiles; i++){
    nn += SifTop[i].nn;  
  }
  fprintf (out,"%ld\n",nn);
  renn = 0;
  for(i = 0; i < nfiles; i++){
    SifTop[i].shift_print_nodes(out,renn);
    renn +=SifTop[i].nn;
  }
  ne = 0;
  for(i = 0; i < nfiles; i++){
    ne += SifTop[i].ne;  
  }
  fprintf (out,"%ld\n",ne);
  
  renn = 0;
  rene = 0;
  for(i = 0; i < nfiles; i++){
    SifTop[i].shift_print_elements(out,renn,rene);
    renn +=SifTop[i].nn;
    rene +=SifTop[i].ne;
  }
  fprintf (out,"\n");
  fprintf (out,"4 %ld\n",nfiles);
  for(i = 0; i < nfiles; i++){
    fprintf (out," %ld ",SifTop[i].nn);
  }
  fprintf (out,"\n");
  long n = 0;
  for(i = 0; i < nfiles; i++){
    for(j = 0; j < SifTop[i].nn; j++){
      fprintf (out,"%ld\n",ltg[n]+1);
      n++;
    }
    fprintf (out,"\n");
  }
  
}


void connectmeshes::run()
{
  //FILE *outtop;
  const char *dot=".";
  char *p;
  char topname[BUFSIZ],topprip[BUFSIZ];
  long i;
  
  
  
  err=read_input_data();
  
  Gtop = new gtopology[nfiles];
  SifTop = new siftop[nfiles];
  Check_unused_prop = 1;
  
  p = strtok(topfile,dot);
  sprintf(topname,"%s",p);
  p = strtok(NULL,dot);
  sprintf(topprip,"%s",p);
  
  //printf("--------------------------\n");
  //printf("READING AND EXPORTING TOPOLOGY\n");
  printf("Reading and exporting topology\n");

  for(i = 0; i < nfiles; i++){
    // input file
    sprintf(d.topf,"%s%ld.%s",topname,i+1,topprip);
    err=read_topology(i);
    if(d.meshfmt == t3d){
      SifTop[i].t3d2sifelformat();
    }
    //printf("Topology from file %s was read and exported\n",d.topf);
    SifTop[i].exporttop (&Gtop[i]);
    Gtop[i].adjacelem (Log);
  } 
  establish_ltg();
  print();  
}



int main(int argc, char *argv[]){
  connectmeshes *connect;
  
  double tws,twe;
  struct tm *local;
  time_t t;
  
  set_prgname(argv[0]);
  printf("-----------------------\n");
  printf("*** MESH CONNECTION ***\n");
  printf("-----------------------\n");
  t = time(NULL);
  local = localtime(&t);
  printf("Program started  %s\n", asctime(local));
  tws = clock();
  connect = new connectmeshes(argc,argv);
  connect->run();
  delete connect;
  twe = clock();
  printf("--------------------------\n");
  printf("Program ended  %s\n", asctime(local));
  printf ("Whole elapsed time of mesh conversion %le s\n",(twe-tws)/(double)CLOCKS_PER_SEC);
  
  
}
