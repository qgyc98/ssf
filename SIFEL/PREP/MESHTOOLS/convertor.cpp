#include "convertor.h"


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
  // section with input files must be detected
  err = xf_setsec(in, bsec_str[begsec_convert]);  
  switch (err)
  {
    case 1:
      print_err("Section with partitioning has not been detected", __FILE__, __LINE__, __func__);
      abort();
    default:
      break;
   }
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
long input_siftop(XFILE *in, descrip *d)
{
  long ret;

  fprintf(stdout, "Reading of mesh topology . . .");
  switch (d->meshfmt)
  {
    case t3d:
      Top->import_t3d(in, d->paral);
      break;
    case sifel:
      ret = Top->read(in, d->paral, d->redgn);
      return(ret);
    default:
      print_err("unknown mesh format is required", __FILE__, __LINE__, __func__);
      return(4);
  }
  return(0);
}



convertor::convertor(int argc, char *argv[])
{
  in = NULL;
  topf = NULL;
  err = 0;
  logname=NULL;
  tmp=NULL;
  err = -1;
  Argc = argc;
  Argv = argv;  


}

convertor::~convertor()
{
  delete Gtm; 
  delete Top;
}



int convertor::read_input_data()
{
  printf("Reading of input data parametrs");
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
  
  
  
  // reading of input files 
  xf_setsec(in, bsec_str[begsec_convert]);
  // reading of line with output file name
  xfscanf(in,"%k%s","output_file_name",&outputfname);
  // conversion type
  xfscanf(in, "%k%m", "conversion_type",&convType_kwdset,&convType);
  // type of processing
  xfscanf(in, "%k%m", "processing",&convprocessing_kwdset,&convprocessing);
  
  switch(convprocessing){
  case seq:{
    //  this is the sequential version
    d.paral=0;
    break;
  }
  case paral:{
    //  this is the parallel version
    d.paral = 1;
    // reading of line with edge number indicator
    xfscanf(in, "%k%ld", "number_of_files", &ndom);
    break;
  }
  }
  xfclose(in);
  printf(" . . . OK\n");
  return(0);
}



int convertor::read_topology()
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
  err = input_siftop(topf,&d);
  xfclose(topf);
  if (err){
    print_err("\nReading of mesh topology failed\n", __FILE__, __LINE__, __func__);
    return(3);
  }
  else{
    printf(" OK\n");
  }
  te = clock();
  printf("Time of reading topology %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  return(0);
  
  
}


int convertor::run()
{
  
  Gtm = new gtopology;
  Top = new siftop;
  Check_unused_prop = 1;
  
  
  err=read_input_data();
  
  switch(convType){
  case T3D2T3D:{
    print_err("Conversion is not yet implemented", __FILE__, __LINE__, __func__);
    break;
  }
  case T3D2SIFEL:{
    printf("-----------------------\n");
    printf("Conversion from T3D to SIFEL format\n");
    if(d.meshfmt == t3d){
      run_t3d2sifel();
    }
    else{
      print_err("Mesh is not in T3D format\n", __FILE__, __LINE__, __func__);
    }
    break;
  }
  case T3D2GID:{
    printf("-----------------------\n");
    printf("Conversion from T3D to GID format\n");
    if(d.meshfmt == t3d){
      run_t3d2gid();
    }
    else{
      print_err("Mesh is not in SIFEL format\n", __FILE__, __LINE__, __func__);
    }
    break;
  }
  case SIFEL2GID:{
    if(d.meshfmt == sifel){
      printf("-----------------------\n");
      printf("Conversion from SIFEL to GID format\n");
      run_sifel2gid();
    }
    else{
      print_err("Mesh is not in SIFEL format\n", __FILE__, __LINE__, __func__);
    }
    break;
  }
  }
  
  return(0);
}


void convertor::run_sifel2gid()
{
  long i,nn,ne;
  FILE *outtop;
  const char *dot=".";
  char *p;
  char pout[BUFSIZ],outname[BUFSIZ],outprip[BUFSIZ],topname[BUFSIZ],topprip[BUFSIZ];
  double ts,te;
  
  
  
  // set tprob for reading of mesh
  switch(convprocessing){
  case seq:{
    printf("Sequential processing\n");
    sprintf(d.topf,"%s",topfile);
    err=read_topology();
    //output file - msh
    p = strtok(outputfname,dot);
    sprintf(pout,"%s.msh",p);
    outtop=fopen(pout,"w");
    printf("Export to GiD format\n");
    ts = clock();
    Top->export_gid_mesh(outtop,1,1);
    te = clock();
    printf("Time of conversion %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
    fclose(outtop);
    break;
  }
  case paral:{
    nn = 1;
    ne = 1;
    //output file
    p = strtok(outputfname,dot);
    sprintf(outname,"%s",p);
    p = strtok(NULL,dot);
    sprintf(outprip,"%s",p);
    p = strtok(topfile,dot);
    sprintf(topname,"%s",p);
    p = strtok(NULL,dot);
    sprintf(topprip,"%s",p);
    printf("Parallel processing\n");
    for(i = 0; i < ndom; i++){
      printf("-----------------------\n");
      printf("Conversion of mesh for subdomain number %ld\n",i+1); 
      // input file
      sprintf(d.topf,"%s%ld.%s",topname,i+1,topprip);
      err=read_topology();
      //output file - msh
      sprintf(pout,"%s%ld.msh",outname,i+1);
      outtop=fopen(pout,"w");
      printf("Export to GiD format\n");
      ts = clock();
      Top->export_gid_mesh(outtop,nn,ne);
      te = clock();
      printf("Time of conversion %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
      fclose(outtop);
      nn+=Top->nn;
      ne+=Top->ne;
    }
    break;
  }
  }
}

void convertor::run_t3d2gid()
{
  long i,nn,ne;
  FILE *outtop;
  const char *dot=".";
  char *p;
  char pout[BUFSIZ],outname[BUFSIZ],outprip[BUFSIZ],topname[BUFSIZ],topprip[BUFSIZ];  
  double ts,te;
    
  switch(convprocessing){
  case seq:{
    printf("Sequential processing\n");
    // input file
    sprintf(d.topf,"%s",topfile);
    err=read_topology();
    Top->t3d2sifelformat();
    //output file - msh
    p = strtok(outputfname,dot);
    sprintf(pout,"%s.msh",p);
    outtop=fopen(pout,"w");
    ts = clock();
    Top->export_gid_mesh(outtop,1,1);
    te = clock();
    printf("Time of conversion %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
    fclose(outtop);
    break;
  }
  case paral:{
    ne = 1;
    nn = 1;
    p = strtok(outputfname,dot);
    sprintf(outname,"%s",p);
    p = strtok(NULL,dot);
    sprintf(outprip,"%s",p);
    p = strtok(topfile,dot);
    sprintf(topname,"%s",p);
    p = strtok(NULL,dot);
    sprintf(topprip,"%s",p);
    printf("Parallel processing\n");
    for(i = 0; i < ndom; i++){
      printf("-----------------------\n");
      printf("Conversion of mesh for subdomain number %ld\n",i+1); 
      // input file
      sprintf(d.topf,"%s.%ld",topfile,i);
      err=read_topology();
      Top->t3d2sifelformat();
      // output file - msh
      sprintf(pout,"%s%ld.msh",outname,i+1);
      outtop=fopen(pout,"w");
      printf("Export to GiD format\n");
      ts = clock();
      Top->export_gid_mesh(outtop,nn,ne);
      te = clock();
      printf("Time of conversion %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
      fclose(outtop);
      nn+=Top->nn;
      ne+=Top->ne;
    }
    break;
  }
  }
}

void convertor::run_t3d2sifel()
{
  long i;
  FILE *outtop;
  const char *dot=".";
  char *p;
  char pout[BUFSIZ],outname[BUFSIZ],outprip[BUFSIZ],topname[BUFSIZ],topprip[BUFSIZ];  
  double ts,te;
  
  switch(convprocessing){
  case seq:{
    printf("Sequential processing\n");
    sprintf(d.topf,"%s",topfile);
    err=read_topology();
    Top->t3d2sifelformat();
    outtop=fopen(outputfname,"w");
    printf("Export to SIFEL format\n");
    ts = clock();
    Top->print(outtop);
    te = clock();
    printf("Time of conversion %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
    fclose(outtop);
    break;
  }
  case paral:{
    p = strtok(outputfname,dot);
    sprintf(outname,"%s",p);
    p = strtok(NULL,dot);
    sprintf(outprip,"%s",p);
    p = strtok(topfile,dot);
    sprintf(topname,"%s",p);
    p = strtok(NULL,dot);
    sprintf(topprip,"%s",p);
    printf("Parallel processing\n");
    for(i = 0; i < ndom; i++){
      printf("-----------------------\n");
      printf("Conversion of mesh for subdomain number %ld\n",i+1); 
      // input file
      sprintf(d.topf,"%s.%ld",topfile,i);
      err=read_topology();
      Top->t3d2sifelformat();
      // output file
      sprintf(pout,"%s%ld.%s",outname,i+1,outprip);
      outtop=fopen(pout,"w");
      printf("Export to SIFEL format\n");
      ts = clock();
      Top->print(outtop);
      te = clock();
      printf("Time of conversion %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
      fclose(outtop);
    }
    break;
  }
  }
}
  


int main (int argc,char *argv[])
{
  
  if (argc < 2){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use :%s input_file_name \n\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  
  convertor *convert;
  double tws,twe;
  struct tm *local;
  time_t t;
  
  set_prgname(argv[0]);
  printf("-----------------------\n");
  printf("*** MESH CONVERSION  ***\n");
  printf("-----------------------\n");
  t = time(NULL);
  local = localtime(&t);
  printf("Program started  %s\n", asctime(local));
  // printf("with this input line parametrs\n",argv[0]);
  //   for(i = 1; i < argc; i++){
  //     printf(" %s",argv[i]);
  //   }
  //   printf("\n");
  
  tws = clock();
  convert = new convertor(argc,argv);
  convert->run();
  delete convert;
  twe = clock();
  printf("--------------------------\n");
  printf("Program ended  %s\n", asctime(local));
  printf ("Whole elapsed time of mesh conversion %le s\n",(twe-tws)/(double)CLOCKS_PER_SEC);
}
