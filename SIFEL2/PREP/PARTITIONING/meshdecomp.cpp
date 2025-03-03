#include "meshdecomp.h"
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
  err = xf_setsec(in, bsec_str[begsec_part]);  
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

  fprintf(stdout, "\nReading of mesh topology . . .");
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


meshdecomp::meshdecomp(int argc, char *argv[])
{
  Argc = argc;
  Argv = argv;  
  in = NULL;
  topf = NULL;
  err = 0;
  logname=NULL;
  tmp=NULL;
}
// descructor
meshdecomp::~meshdecomp()
{
  //delete Gtm;
  delete Top;
  
}

int meshdecomp::read_input_data()
{
  printf("--------------------------\n");
  printf("DATA READING\n");
  printf("\nReading of input data parametrs");
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
  xfscanf(in, "%k%s","topology_file",d.topf);
  // reading of line with topology file format indicator,
  xfscanf(in, "%k%m", "mesh_format", &meshform_kwdset, &d.meshfmt);
  // reading of line with edge number indicator
  xfscanf(in, "%k%ld", "edge_numbering", &d.redgn);
  //  this is the sequential version
  d.paral=0;
  
  // reading of input files 
  xf_setsec(in, bsec_str[begsec_part]);
  // type of processing
  xfscanf(in, "%k%m", "processing",&processing_kwdset,&processing);
  
  
  switch(processing){
  case preprocessing:{
    // partitioning type
    xfscanf(in, "%k%m", "partType",&partType_kwdset,&partType);
    // weightred graph
    xfscanf(in, "%k%m", "weight_graph",&graphWeight_kwdset,&graphWeight);
    if(graphWeight != noweight){
      xfscanf(in, "%k%m", "weight_option",&partOption_kwdset,&weightOpt);
      if(weightOpt == user){
	xfscanf(in, "%k%a","file_with_weights",weightf);
      }
    }
    xfscanf(in, "%k%s","graph_output_file",&graphf);
    break;
  }
  case postprocessing:{
    //  type of mesh description
    xfscanf (in,"%k%m","mesh_description",&meshdescription_kwdset,(int*)&md);
    // partitioning type
    xfscanf(in, "%k%m", "partType",&partType_kwdset,&partType);
    // number of partition
    xfscanf(in, "%k%ld", "number_of_partitions", &nparts);
    // reading file with partitioning
    xfscanf(in, "%k%s","partitioning_file",&partf);
    // reading of line with output file name
    xfscanf(in,"%k%s","output_file_name",&outputfname);
    break;
  }
  case all:{
    //  type of mesh description
    xfscanf (in,"%k%m","mesh_description",&meshdescription_kwdset,(int*)&md);
    // // partitioning type
    //     xfscanf(in, "%k%m", "partType",&partType_kwdset,&partType);
    partType=METIS;
    // number of partition
    xfscanf(in, "%k%ld", "number_of_partitions", &nparts);
    // weightred graph
    xfscanf(in, "%k%m", "weight_graph",&graphWeight_kwdset,&graphWeight);
    if(graphWeight != noweight){
      xfscanf(in, "%k%m", "weight_option",&partOption_kwdset,&weightOpt);
      if(weightOpt == user){
	xfscanf(in, "%k%a","file_with_weights",weightf);
      }
    }
    // reading of line with output file name
    xfscanf(in,"%k%s","output_file_name",&outputfname);
    switch(partType){
    case METIS:{
      read_metis_option();
      break;
    }
    case CHACO:{
      
      break;
    }
    case JOSTLE:{
      
      break;
    }
    case PARTY:{
      
      break;
    }
    case SCOTCH:{
      
      break;
    }
    }
    break;
  }
  case aggregates:{
    partType = METIS;
    // number of partition
    xfscanf(in, "%k%ld", "number_of_aggregates", &nparts);
    graphWeight = noweight;
    // weightred graph
    //  xfscanf(in, "%k%m", "weight_graph",&graphWeight_kwdset,&graphWeight);
    //     if(graphWeight != noweight){
    //       xfscanf(in, "%k%m", "weight_option",&partOption_kwdset,&weightOpt);
    //       if(weightOpt == user){
    // 	xfscanf(in, "%k%a","file_with_weights",weightf);
    //       }
    //     }
    // reading of line with output file name
    xfscanf(in,"%k%s","output_file_name",&outputfname);
    read_metis_option();
    break;
  }
  }
  
  
  xfclose(in);
  printf(" . . . OK");
  return(0);
}



int meshdecomp::read_topology()
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
    printf(" OK");
  }
  te = clock();
  printf("\nTime of reading topology %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  return(0);

  // if T3D topology -> check and renumber
  if(d.meshfmt == t3d){
    Top->t3d2sifelformat();
  }
  
}

void meshdecomp::read_metis_option()
{
  xfscanf(in, "%k%m", "partitioning_technique",&partTech_kwdset,&parttech);
  switch(parttech){
  case recursive:{
    //xfscanf(in, "%k%m", "partitioning_option",&partOption_kwdset,&opt);
    options = new int[5];
    opt = defaul;
    if(opt == user){
      options[0]=1;
    // future time
    }
    else{
      options[0]=0;
    }
    break;
  }
  case kway:{
    //xfscanf(in, "%k%m", "partitioning_option",&partOption_kwdset,&opt);
    options = new int[5];
    opt = defaul;
    if(opt == user){
      options[0]=1;
      // future time
    }
    else{
      options[0]=0;
    }  
    break;
  } 
  case vkway:{
    //xfscanf(in, "%k%m", "partitioning_option",&partOption_kwdset,&opt);
    opt = defaul;
    options = new int[5];
    if(opt == user){
      options[0]=1;
      // future time
    }
    else{
    options[0]=0;
    } 
    break;
  }
  case mcrekursive:{
    xfscanf(in, "%k%d", "number_of_constraints",&partOption_kwdset,&ncon);
    //xfscanf(in, "%k%m", "partitioning_option",&partOption_kwdset,&opt);
    options = new int[5];
    opt = defaul;
    if(opt == user){
      options[0]=1;
      // future time
    }
    else{
    options[0]=0;
    } 
    break;
  }
  case mckway:{
    xfscanf(in, "%k%d", "number_of_constraints",&partOption_kwdset,&ncon);
    xfscanf(in, "%k%m", "file_ncon_tolerance",&ubecfile);
    //xfscanf(in, "%k%m", "partitioning_option",&partOption_kwdset,&opt);
    opt = defaul;
    options = new int[5];
    if(opt == user){
      options[0]=1;
      // future time
    }
    else{
      options[0]=0;
    } 
    break;
  }
  case wpartrekursive:{
    xfscanf(in, "%k%m", "file_partition_weights",&tpwgts);
    //xfscanf(in, "%k%m", "partitioning_option",&partOption_kwdset,&opt);
    opt = defaul;
    options = new int[5];
    if(opt == user){
      options[0]=1;
      // future time
    }
    else{
      options[0]=0;
    }    
    break;
  }
  case wpartkway:{
    xfscanf(in, "%k%m", "file_partition_weights",&tpwgts);
    //xfscanf(in, "%k%m", "partitioning_option",&partOption_kwdset,&opt);
    opt = defaul;
    options = new int[5];
    if(opt == user){
      options[0]=1;
      // future time
    }
    else{
      options[0]=0;
    } 
    break;
  }
  case wpartvkway:{
    xfscanf(in, "%k%m", "file_partition_weights",&tpwgts);
    //xfscanf(in, "%k%m", "partitioning_option",&partOption_kwdset,&opt);
    opt = defaul;
    options = new int[5];
    if(opt == user){
      options[0]=1;
      // future time
    }
    else{
      options[0]=0;
    } 
    break;
  }
  }  
  
  numflag = 0;
  switch(graphWeight){
  case noweight:{
    wgtflag = 0;
    break;
  }
  case nodalweight:{
    wgtflag = 2;
    break;
  }
  case edgeweight:{
    wgtflag = 1;
    break;
  }
  case nodealedgeweight:{
    wgtflag = 3;
    break;
  }
  }
  

  
}

void meshdecomp::read_chaco_option()
{

}


int meshdecomp::run_metis()
{
  double ts,te;
  long i;
  switch(processing){
  case preprocessing:{
    printf("--------------------------\n");
    printf("TOPOLOGY EXPORT\n");
    Top->exporttop (Gtm);
    printf("--------------------------\n");
    creation_dual_graph();
    print_metis_graph();
    break;
  }
  case postprocessing:{
    read_partitioning();
    part_postprocessing();
    break;
  }
  case all:{
    printf("--------------------------\n");
    printf("TOPOLOGY EXPORT\n");
    Top->exporttop (Gtm);
    printf("--------------------------\n");

    creation_dual_graph();
    switch(parttech){
    case mckway:{
      ubec = new float[ncon];
      read_metis_ubec();
      break;
    }
    case wpartrekursive:{
      tpwgts = new float[nparts];
      read_metis_tpwgts();
      break;
    }
    case wpartkway:{
      tpwgts = new float[nparts];
      read_metis_tpwgts();
      break;
    }
    case wpartvkway:{
      tpwgts = new float[nparts];
      read_metis_tpwgts();
      break;
    }
    }
    
    printf("--------------------------\n");
    printf("METIS PARTITIONING\n");
    printf("Graph Information\n");
    printf("#Vertices: %d, #Edges: %d, #Parts: %d\n",nvtxs, nedges/2, nparts);
    
    metis_part = new idxtype[nvtxs];
    metis_nvtxs = nvtxs;
    metis_xadj = new idxtype[metis_nvtxs+1];
    if(vwgt != NULL)  metis_vwgt = new idxtype[metis_nvtxs];
    for(i = 0; i < metis_nvtxs; i++){
      metis_xadj[i]=xadj[i];
      if(vwgt != NULL) metis_vwgt[i]=vwgt[i];
    }
    metis_xadj[metis_nvtxs]=xadj[nvtxs];
    metis_adjncy = new idxtype[nedges]; 
    if(adjwgt != NULL) metis_adjwgt = new idxtype[nedges]; 
    for(i = 0; i < nedges; i++){
      metis_adjncy[i] = adjncy[i];
      if(adjwgt != NULL)  metis_adjwgt[i] = adjwgt[i];
    }
    delete []xadj;
    delete []adjncy;
    delete []vwgt;
    delete []adjwgt;

    // due to metis manual
    if(nparts < 8 && parttech != recursive){
      options[0] = 0;
      parttech = recursive;
      printf("Because nparts < 8 Kway changed to resursive bisection due to recommendation in metis manual\n");
    }
    ts = clock();
    // promene definovane na miste je treba naplnit pri cteni promenich ze sekce partitioning
    switch(parttech){
    case recursive:{
      METIS_PartGraphRecursive(&metis_nvtxs, metis_xadj, metis_adjncy, metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, options, &edgecut, metis_part);
      break;
    }
    case kway:{
      METIS_PartGraphKway(&metis_nvtxs, metis_xadj, metis_adjncy, metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, options, &edgecut, metis_part);
      //printf("edge cut %ld\n",edgecut);
      break;
    } 
    case vkway:{
      METIS_PartGraphVKway(&metis_nvtxs, metis_xadj, metis_adjncy, metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, options, &volume, metis_part);
      printf("volume %d\n",volume);
      break;
    }
    case mcrekursive:{
      METIS_mCPartGraphRecursive(&metis_nvtxs,&ncon, metis_xadj, metis_adjncy, metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, options, &edgecut, metis_part);
      break;
    }
    case mckway:{
      METIS_mCPartGraphKway(&metis_nvtxs,&ncon, metis_xadj, metis_adjncy, metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, ubec,options, &edgecut, metis_part);
      break;
    }
    case wpartrekursive:{
      METIS_WPartGraphRecursive(&metis_nvtxs, metis_xadj, metis_adjncy, metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, tpwgts,options, &edgecut, metis_part);
      break;
    }
    case wpartkway:{
      METIS_WPartGraphKway(&metis_nvtxs, metis_xadj, metis_adjncy,metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, tpwgts,options, &edgecut, metis_part);
      break;
    }
    case wpartvkway:{
      METIS_WPartGraphVKway(&metis_nvtxs, metis_xadj, metis_adjncy, metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, tpwgts,options, &volume, metis_part);
      break;
    }
    }
    te = clock();
    printf("Time of metis partitoning %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
    
    delete []metis_xadj;
    delete []metis_adjncy;
    delete []metis_vwgt;
    delete []metis_adjwgt;
    
    part = new int[nvtxs];
    for(i = 0; i < nvtxs; i++){
      part[i] = metis_part[i];
    }
    delete []metis_part;
    part_postprocessing();
    
    break;
  }
  }
  return 1;
}

void meshdecomp::read_metis_ubec()
{
  FILE *in;
  long i;
  in = fopen(ubecfile,"r");
  if (in == NULL){
    fprintf (stderr,"\n File with ubec vector specification has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
  }
  for(i = 0; i < ncon; i++){
    fscanf(in,"%f",&ubec[i]);
  }
  fclose(in);
  
}

void meshdecomp::read_metis_tpwgts()
{
  FILE *in;
  long i;
  in = fopen(tpwgtsfile,"r");
  if (in == NULL){
    fprintf (stderr,"\n File with tpwgts vector sepecification has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
  }
  for(i = 0; i < nparts; i++){
    fscanf(in,"%f",&tpwgts[i]);
  }
  fclose(in);
 
}


void meshdecomp::print_metis_graph()
{
  FILE *out;
  long i,j;
  long adr1,adr2;
  char outt[1050]; 
  
  sprintf(outt,"%s",graphf);
  out = fopen(outt,"w");
  if (out == NULL){
    fprintf (stderr,"\n File for metis graph output has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
  }

  switch(weightOpt){
  case noweight:{
    fprintf(out,"%d %d\n",nvtxs,nedges/2);
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d",adjncy[j]+1);
      }
      fprintf(out,"\n");
    }
    break;
  }
  case nodalweight:{
    fprintf(out,"%d %d 10 1\n",nvtxs,nedges/2);
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      fprintf(out," %d",vwgt[i]);
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d",adjncy[j]+1);
      }
      fprintf(out,"\n");
    }

    break;
  }
  case edgeweight:{
    fprintf(out,"%d %d 1\n",nvtxs,nedges/2);
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d   %d",adjncy[j]+1,adjwgt[j]);
      }
      fprintf(out,"\n");
    }
    break;
  }
  case nodealedgeweight:{
    fprintf(out,"%d %d 1\n",nvtxs,nedges/2);
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      fprintf(out," %d",vwgt[i]);
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d   %d",adjncy[j]+1,adjwgt[j]);
      }
      fprintf(out,"\n");
    }
    break;
  }
  case multiconstrain:{
    fprintf(out,"%d %d 10 %d\n",nvtxs,nedges/2,ncon);
    for(i = 0; i < nvtxs; i++){
      for(j = 0; j < ncon; j++){
	fprintf(out,"   %d",vwgt[i+j]);
      }
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d",adjncy[j]+1);
      }
      fprintf(out,"\n");
    }
    break;
  }
  }
  fclose(out);
  
}


void meshdecomp::print_chaco_graph()
{
  FILE *out;
  long i,j;
  long adr1,adr2;
  char outt[1050]; 
  
  sprintf(outt,"%s.graph",graphf);
  out = fopen(outt,"w");
  if (out == NULL){
    fprintf (stderr,"\n File for chaco graph output has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
  }
  switch(weightOpt){
  case noweight:{
    fprintf(out,"%d %d\n",nvtxs,nedges/2);
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d",adjncy[j]+1);
      }
      fprintf(out,"\n");
    }
    break;
  }
  case nodalweight:{
    fprintf(out,"%d %d 101\n",nvtxs,nedges/2);
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      fprintf(out,"%ld  %d",i+1,vwgt[i]);
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d",adjncy[j]+1);
      }
      fprintf(out,"\n");
    }

    break;
  }
  case edgeweight:{
    fprintf(out,"%d %d  11\n",nvtxs,nedges/2);
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      fprintf(out,"%ld   ",i+1);
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d   %d",adjncy[j]+1,adjwgt[j]);
      }
      fprintf(out,"\n");
    }
    break;
  }
  case nodealedgeweight:{
    fprintf(out,"%d %d 111\n",nvtxs,nedges/2);
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      fprintf(out,"%ld   %d",i+1,vwgt[i]);
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d   %d",adjncy[j]+1,adjwgt[j]);
      }
      fprintf(out,"\n");
    }
    break;
  }
  }
  fclose(out);
  
  sprintf(outt,"%s.coords",graphf);

  out = fopen(outt,"w");
  if (out == NULL){
    fprintf (stderr,"\n File for chaco coords output has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
  }
  for(i = 0; i < nvtxs; i++){
    fprintf(out,"  %lf   %lf %lf\n",x[i],y[i],z[i]);
  }
  fclose(out);
}
void meshdecomp::print_jostle_graph()
{
  FILE *out;
  long i,j;
  long adr1,adr2;
  char outt[1050]; 
  
  sprintf(outt,"%s.grf",graphf);
  out = fopen(outt,"w");
  if (out == NULL){
    fprintf (stderr,"\n File for jostle graph output has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
  }
  switch(weightOpt){
  case noweight:{
    fprintf(out,"%d %d\n",nvtxs,nedges/2);
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d",adjncy[j]+1);
      }
      fprintf(out,"\n");
    }
    break;
  }
  case nodalweight:{
    fprintf(out,"%d %d 10\n",nvtxs,nedges/2);
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      fprintf(out," %d",vwgt[i]);
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d",adjncy[j]+1);
      }
      fprintf(out,"\n");
    }

    break;
  }
  case edgeweight:{
    fprintf(out,"%d %d 1\n",nvtxs,nedges/2);
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d   %d",adjncy[j]+1,adjwgt[j]);
      }
      fprintf(out,"\n");
    }
    break;
  }
  case nodealedgeweight:{
    fprintf(out,"%d %d 11\n",nvtxs,nedges/2);
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      fprintf(out," %d",vwgt[i]);
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d   %d",adjncy[j]+1,adjwgt[j]);
      }
      fprintf(out,"\n");
    }
    break;
  }
  }
  fclose(out);
}



void meshdecomp::print_scotch_graph()
{
  FILE *out;
  long i,j;
  long adr1,adr2;
  char outt[1050]; 
  
  sprintf(outt,"%s.grf",graphf);
  out = fopen(outt,"w");
  if (out == NULL){
    fprintf (stderr,"\n File for scotch graph output has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
  }
  fprintf(out,"0\n");
  fprintf(out,"%d %d\n",nvtxs,nedges);
  switch(weightOpt){
  case noweight:{
    fprintf(out,"0  000\n");
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      fprintf(out,"  %ld  ",adr2-adr1);
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d",adjncy[j]+1);
      }
      fprintf(out,"\n");
    }
    break;
  }
  case nodalweight:{
    fprintf(out,"0 101\n");
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      fprintf(out,"   %d",vwgt[i]);
      fprintf(out,"  %ld  ",adr2-adr1);
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d",adjncy[j]+1);
      }
      fprintf(out,"\n");
    }

    break;
  }
  case edgeweight:{
    fprintf(out,"0  011\n");
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      fprintf(out,"  %ld  ",adr2-adr1);
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d   %d",adjncy[j]+1,adjwgt[j]);
      }
      fprintf(out,"\n");
    }
    break;
  }
  case nodealedgeweight:{
    fprintf(out,"0 111\n");
    for(i = 0; i < nvtxs; i++){
      adr1 = xadj[i];
      adr2 = xadj[i+1];
      fprintf(out,"%d",vwgt[i]);
      fprintf(out,"  %ld  ",adr2-adr1);
      for(j = adr1; j < xadj[i+1]; j++){
	fprintf(out,"  %d   %d",adjncy[j]+1,adjwgt[j]);
      }
      fprintf(out,"\n");
    }
    break;
  }
  }
  fclose(out);
  
  sprintf(outt,"%s.xyz",graphf);
  out = fopen(outt,"w");
  if (out == NULL){
    fprintf (stderr,"\n File for scotch coords output has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
  }
  fprintf(out,"3\n%d\n",nvtxs);
  for(i = 0; i < nvtxs; i++){
    fprintf(out,"%ld  %lf   %lf\n",i,x[i],y[i]);
  }
  fclose(out);
}

void meshdecomp::compute_centre_elements()
{
  long i,j,k;
  x = new double[nvtxs];
  y = new double[nvtxs];
  z = new double[nvtxs];
  
  
  for(i = 0; i < Top->ne; i++){
    x[i] = 0.0;
    y[i] = 0.0;
    z[i] = 0.0;
    for(j = 0; j < Top->elements[i].nne; j++){ 
      k=Top->elements[i].nodes[j];
      x[i]+=Top->nodes[k].x;
      y[i]+=Top->nodes[k].y;
      z[i]+=Top->nodes[k].z;
    }
    x[i]/=Top->elements[i].nne;
    y[i]/=Top->elements[i].nne;
    z[i]/=Top->elements[i].nne;
  }
  
}

void meshdecomp::read_partitioning()
{
  FILE *in;
  long i;
  nvtxs = Top->ne;
  in = fopen(partf,"r");
  if (in == NULL){
    fprintf (stderr,"\n File with partitoning has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
  }
  part = new int[nvtxs];
  for(i = 0; i < nvtxs; i++){
    fscanf(in,"%d",&part[i]);
  }
  fclose(in);
}

void  meshdecomp::part_postprocessing()
{
  long i,j,n,ii,jj,nnd,k,l;
  long nbn;
  long *lned,*nod,*multip,*nnsd;
  //  long **part_nodes;
  char outt[2*BUFSIZ+50]; 
  FILE *out;
  
  printf("--------------------------\n");
  printf("POSTPROCESSING OF PARTITIONING\n");
  
  //  nodal multiplicity 
  multip = new long [Top->nn];
  for(i = 0; i < Top->nn; i++){
    multip[i]=0;
  }
  
  nod = new long [Top->nn];
  nnsd = new long[nparts];
  lned = new long [nparts];
  
  for (i = 0; i < nparts; i++){
    nnsd[i] = 0;
    lned[i] = 0;
    
    for (j = 0; j < Top->nn; j++){
      nod[j]=0;
    }
    
    // number of elements on the i-th subdomain
    for (j = 0;j < Top->ne; j++){
      if (part[j] == i){
	lned[i]++;
	for(k = 0; k < Top->elements[i].nne; k++){
	  jj=Top->elements[j].nodes[k];
	  nod[jj]++;
	}
      }
    }
    
    //     for (j = 0; j < Top->nn; j++){
    //       printf("nod[%ld]=%ld\n",j+1,nod[j]);
    //     }
    
    // number of nodes on the i-th subdomain
    nnd=0;
    for (j = 0; j < Top->nn; j++){
      if (nod[j]>0){
	nnd++;
	nnsd[i]++;
	multip[j]++;
      }
    }
  }
  
  // for (j = 0; j < Top->nn; j++){
  //     printf("multip[%ld]=%ld\n",j+1,multip[j]);
  //   }
  

  
  printf("Number of elements and number of nodes on particular subdomains\n"); 
  
  for (i=0;i<nparts;i++){
    printf("Sudomain # %ld : # elements  %ld  #nodes %ld\n",i,lned[i],nnsd[i]);
  }
  
  // check of sum of elements
  n = 0;
  for(i = 0; i < nparts; i++){ 
    n+=lned[i];
  }
  if (n != Top->ne)
    printf("Sum of elements on subdomains is not equal to the number of elements in original mesh!\n");
  else
    printf("Sum of elements on subdomains is equal to the number of elements in original mesh\n"); 


  // mesh description
  switch(md){
  case all_nodes:{
    nbn=1;
    for (i = 0;i < Top->nn; i++){
      if (multip[i] > 1){
	multip[i]=i+1; 
	nbn++;
      }
      else{
	multip[i]=i+1; 
      }
    }
    break;
  }
    // boundary nodes
  case bound_nodes:{
    nbn=1;
    for (i = 0; i < Top->nn; i++){
      if (multip[i] > 1){
	multip[i] = nbn;
	nbn++;
      }
      else{
	multip[i] = 0;
      }
    }
    break;
  }
  case neg_bound_nodes:{
    nbn=1;
    for (i = 0;i < Top->nn; i++){
      if (multip[i] > 1){
	multip[i] = (i+1)*(-1); 
	nbn++;
      }
      else{
	multip[i] = i+1; 
      }
    }
    break;
  }
  case  metis:{
    break;
  }
  case glob_glued:{
     
    break;
  }
  default:
    fprintf(stderr, "Unknown type of mesh description (md=%d)\n", md);
    abort();
  }
  printf ("Number of boundary nodes  %ld\n",nbn);
  
  
  // for (j = 0; j < Top->nn; j++){
  //       printf("multip[%ld]=%ld\n",j+1,multip[j]);
  //     }
  const char *dot=".";
  char *p;
  char name[BUFSIZ],prip[BUFSIZ];
  //output file
  p = strtok(outputfname,dot);
  sprintf(name,"%s",p);
  p = strtok(NULL,dot);
  sprintf(prip,"%s",p);
  for (i = 0; i < nparts; i++){
    sprintf(outt,"%s%ld.%s",name,i+1,prip);
    printf("Writting of topology of subdomain # %ld into file %s\n",i,outt);
    out=fopen(outt,"wt");
    if (out == NULL){
      fprintf (stderr,"\n File for output topology has not been specified.");
      fprintf (stderr,"\n Try it again!\n\n");
    }
    for (j = 0; j < Top->nn; j++){
      nod[j]=0;
    }
    
    for (j = 0;j < Top->ne; j++){
      if (part[j] == i){
	for(k = 0; k < Top->elements[i].nne; k++){
	  jj=Top->elements[j].nodes[k];
	  nod[jj]++;
	}
      }
    }
    //     for (j = 0; j < Top->nn; j++){
    //       printf("nod[%ld]=%ld\n",j+1,nod[j]);
    //     }
    

    
//    fprintf (out,"# Parallel topology\n");
//    fprintf (out,"%ld # number of nodes\n",nnsd[i]);
    fprintf (out,"%ld\n",nnsd[i]);
   
    // printing nodes
    k = 1;
    for (ii = 0; ii < Top->nn; ii++){
      if (nod[ii] > 0){
	fprintf (out,"%ld %15.10le %15.10le %15.10le %ld",k, Top->nodes[ii].x, Top->nodes[ii].y, Top->nodes[ii].z, Top->nodes[ii].nprop);
	for(j = 0; j < Top->nodes[ii].nprop; j++){
	  fprintf(out, " %d %ld", Top->nodes[ii].entid[j], Top->nodes[ii].prop[j]);
	}
	fprintf(out, "\n");
	nod[ii]=k;
	k++;
      }
    }
//     for (j = 0; j < Top->nn; j++){
//       printf("nod[%ld]=%ld\n",j+1,nod[j]);
//     }
    // printing elements 
//    fprintf (out,"\n%ld #number of elements\n", lned[i]);
    fprintf (out,"\n%ld\n", lned[i]);
    l = 1;
    for (ii = 0; ii < Top->ne; ii++){
      if (part[ii]==i){
	fprintf (out,"%ld %d",l, Top->elements[ii].type);
	for (j = 0; j < Top->elements[ii].nne; j++){
	  fprintf (out," %ld",nod[Top->elements[ii].nodes[j]]);
	}
	fprintf (out," %ld",Top->elements[ii].prop);
	if (Top->elements[ii].propedg){
	  fprintf(out, "  ");
	  for (j=0; j < Top->elements[ii].ned; j++)
	    fprintf (out," %ld",Top->elements[ii].propedg[j]);
	}
	if (Top->elements[ii].propsurf){
	  fprintf(out, "  ");
	  for (j=0; j < Top->elements[ii].nsurf; j++)
	    fprintf (out," %ld",Top->elements[ii].propsurf[j]);
	}
	fprintf(out, "\n");
	l++;
      }
    }
    // printings ltg
    for (j = 0; j < Top->nn; j++){
      if (nod[j] > 0){
	fprintf (out,"%5ld   %5ld\n",nod[j],multip[j]);
      }
    }
    
    
    // printing edges
    if (Top->edges)
      Top->edges->print(out, nod);
    // printing surfaces
    if (Top->surfaces)
      Top->surfaces->print(out, nod);

    fclose (out);
  }
  
  


  // delete []nod;
  //delete []multip;
  //delete []lned;
  
}

void meshdecomp::aggregate_postprocessing()
{
  printf("--------------------------\n");
  printf("POSTPROCESSING OF PARTITIONING\n");
  
  long i,j;
  long *nnag;  
  FILE *topfile;
  
  nnag = new long[nparts];
  for(i = 0; i < nparts; i++){
    nnag[i]=0; 
  }
  for(i = 0; i < Top->nn; i++){
    j = part[i];
    nnag[j]++;
  }
  topfile=fopen(outputfname,"w");
  Top->print (topfile);
  
  fprintf(topfile,"\n");
  fprintf(topfile,"11 %d\n",nparts);
  for(i = 0; i < nparts; i++){
    fprintf(topfile,"%ld\n",nnag[i]);
  }
  fprintf(topfile,"\n");
  for(i = 0; i < Top->nn; i++){
    fprintf(topfile,"%d\n",part[i]);
  }

  fclose(topfile);
  
  


}

void meshdecomp::creation_dual_graph()
{
  long i,j,k;
  double ts,te;
  ts = clock();
  printf("--------------------------\n");
  printf("DUAL GRAPH CREATION\n");
  Gtm->adjacelem (Log);
  
  nvtxs = Top->ne;

  if (nvtxs <= 0) {
    printf("Empty  Nothing to do.\n");
    exit(0);
  }

  xadj = new int[nvtxs+1];
  
  xadj[0] = 0;
  
  j = 0;
  for(i = 1; i < nvtxs; i++){
    j+=Gtm->nadjelel[i-1]-1;
    xadj[i] = j;
  }
  j+=Gtm->nadjelel[nvtxs-1]-1;
  xadj[nvtxs]=j;
  nedges = xadj[nvtxs];
  adjncy = new int[nedges];
  k = 0;
  for(i = 0; i < nvtxs; i++){
    for(j = 0; j < Gtm->nadjelel[i]; j++){
      if(Gtm->adjelel[i][j] != i){
	adjncy[k] = Gtm->adjelel[i][j];
	k++;
      }
    }
  }
  delete Gtm;
  
  if(graphWeight != noweight){
    if(weightOpt == user){
      //read_metis_weight();
    }
    else{
      //compute_metis_weight();
    }
  }
  else{
    // vertex weight
    vwgt = NULL;
    //edge weight
    adjwgt = NULL;
  }
  te = clock();
  printf("Time of creation of dual graph %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  
  
}

void meshdecomp::creation_nodal_graph()
{

  long i,j,k;
  double ts,te;
  ts = clock();
  printf("NODAL GRAPH CREATION\n");
  
  
//   Gtm->adjacnodes (Log);
     
//   nvtxs = Top->nn;

//   if (nvtxs <= 0) {
//     printf("Empty  Nothing to do.\n");
//     exit(0);
//   }

//   xadj = new int[nvtxs+1];
  
//   j = 0;
//   for(i = 1; i < nvtxs; i++){
//     j+=Gtm->nadjnodnod[i-1]-1;
//     xadj[i] = j;
//   }
//   j+=Gtm->nadjnodnod[nvtxs-1]-1;
//   xadj[nvtxs]=j;
//   nedges = xadj[nvtxs];
//   adjncy = new int[nedges];
//   k = 0;
//   for(i = 0; i < nvtxs; i++){
//     for(j = 0; j < Gtm->nadjnodnod[i]; j++){
//       if(Gtm->adjnodnod[i][j] != i){
// 	adjncy[k] = Gtm->adjnodnod[i][j];
// 	k++;
//       }
//     }
//   }

  
  Gtm->adjacnodes_edge (Log);
  
    nvtxs = Gtm->nn;
  
    if (nvtxs <= 0) {
      printf("Empty  Nothing to do.\n");
      exit(0);
    }
  
    xadj = new int[nvtxs+1];
  
  
    xadj[0] = 0;
    j = 0;
    for(i = 1; i < nvtxs; i++){
      j+=Gtm->nadjacnodesedge[i];
      xadj[i] = j;
    }
    j+=Gtm->nadjacnodesedge[nvtxs-1];
    xadj[nvtxs]=j;
    nedges = xadj[nvtxs];
    adjncy = new int[nedges];
    k = 0;
    for(i = 0; i < nvtxs; i++){
      for(j = 0; j < Gtm->nadjacnodesedge[i]; j++){
        adjncy[k] = Gtm->adjacnodesedge[i][j];
        k++;
      }
    }
    
    delete Gtm;
  
  if(graphWeight != noweight){
    if(weightOpt == user){
      //read_metis_weight();
    }
    else{
      //compute_metis_weight();
    }
  }
  else{
    // vertex weight
    vwgt = NULL;
    //edge weight
    adjwgt = NULL;
  }
  te = clock();
  printf("Time of creation of nodal graph %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  
  
}

int meshdecomp::run_chaco()
{
  
  switch(processing){
  case preprocessing:{
    printf("Chaco preprocessing\n");
    printf("--------------------------\n");
    printf("TOPOLOGY EXPORT\n");
    Top->exporttop (Gtm);
    printf("--------------------------\n");

    creation_dual_graph();
    compute_centre_elements();
    print_chaco_graph();
    break;
  }
  case postprocessing:{
    read_partitioning();
    part_postprocessing();
    break;
  }
  case all:{
    printf("Chaco - all - Sorry not supported option in this time\n");
    break;
  }
  }
    
  return 1;
}

void meshdecomp::run_jostle()
{
  switch(processing){
  case preprocessing:{
    printf("--------------------------\n");
    printf("TOPOLOGY EXPORT\n");
    Top->exporttop (Gtm);
    printf("--------------------------\n");

    creation_dual_graph();
    print_jostle_graph();
    break;
  }
  case postprocessing:{
    read_partitioning();
    part_postprocessing();
    break;
  }
  case all:{
    printf("Jostle - all - Sorry not supported option in this time\n");
    break;
  }
  }

}

void meshdecomp::run_scotch()
{
  switch(processing){
  case preprocessing:{
    printf("--------------------------\n");
    printf("TOPOLOGY EXPORT\n");
    Top->exporttop (Gtm);
    
    printf("--------------------------\n");
    creation_dual_graph();
    compute_centre_elements();
    print_scotch_graph();
    break;
  }
  case postprocessing:{
    read_partitioning();
    part_postprocessing();
    break;
  }
  case all:{
    printf("Scotch - all - Sorry not supported option in this time\n");
    break;
  }
  }

}

int meshdecomp::run_party()
{

  // double ts,te;
//   int  cutsize, redl=0,  Output=1, Times=2, rec=1, n=0, *edge_p=NULL, *edge=NULL, *edge_w=NULL, *part;
//   char        	global[100]="def", local[100]="hs",redm[100]="lam", redo[100]="w3";
//   float *xx,*yy,*zz;
  
  
  switch(processing){
  case preprocessing:{
    printf("--------------------------\n");
    printf("TOPOLOGY EXPORT\n");
    Top->exporttop (Gtm);
    printf("--------------------------\n");

    creation_dual_graph();
    compute_centre_elements();
    print_chaco_graph();
    break;
  }
  case postprocessing:{
    read_partitioning();
    part_postprocessing();
    break;
  }
  case all:{
    printf("Party - all - Sorry not supported option in this time\n");
    printf("--------------------------\n");
    printf("TOPOLOGY EXPORT\n");
    Top->exporttop (Gtm);
    
    printf("--------------------------\n");

    creation_dual_graph();
    compute_centre_elements();
    //int err;
    //err=party_lib(nvtxs,vwgt,xx,yy,zz,xadj,adjncy,adjwgt,nparts,part,&cutsize,redl,redm,redo,global,local,rec,Output-1);
    break;
  }
  }
    
  return 1;
}

void meshdecomp::run_aggregates()
{
  long i;
  double ts,te;
  
  printf("--------------------------\n");
  printf("TOPOLOGY EXPORT\n");
  Top->exporttop (Gtm);
  
  printf("--------------------------\n");

  creation_nodal_graph();

  switch(parttech){
  case mckway:{
    ubec = new float[ncon];
    read_metis_ubec();
    break;
  }
  case wpartrekursive:{
    tpwgts = new float[nparts];
    read_metis_tpwgts();
    break;
  }
  case wpartkway:{
    tpwgts = new float[nparts];
    read_metis_tpwgts();
    break;
  }
  case wpartvkway:{
    tpwgts = new float[nparts];
    read_metis_tpwgts();
    break;
  }
  }
  
  printf("--------------------------\n");
  printf("METIS PARTITIONING\n");
  printf("Graph Information\n");
  printf("#Vertices: %d, #Edges: %d, #Parts: %d\n",nvtxs, nedges/2, nparts);
  
  metis_part = new idxtype[nvtxs];
  metis_nvtxs = nvtxs;
  metis_xadj = new idxtype[metis_nvtxs+1];
  if(vwgt != NULL){
    metis_vwgt = new idxtype[metis_nvtxs];
  }
  for(i = 0; i < metis_nvtxs; i++){
    metis_xadj[i]=xadj[i];
    if(vwgt != NULL){
      metis_vwgt[i]=vwgt[i];
    }
  }
  metis_xadj[metis_nvtxs]=xadj[nvtxs];
  metis_adjncy = new idxtype[nedges]; 
  if(adjwgt != NULL){
    metis_adjwgt = new idxtype[nedges]; 
  }
  for(i = 0; i < nedges; i++){
    metis_adjncy[i] = adjncy[i];
    if(adjwgt != NULL){
      metis_adjwgt[i] = adjwgt[i];
    }
  }
  delete []xadj;
  delete []adjncy;
  delete []vwgt;
  delete []adjwgt;
  
  // due to metis manual
  if(nparts < 8 && parttech != recursive){
    options[0] = 0;
    parttech = recursive;
    printf("Because nparts < 8 Kway changed to resursive bisection due to recommendation in metis manual\n");
  }
  ts = clock();
  // promene definovane na miste je treba naplnit pri cteni promenich ze sekce partitioning
  switch(parttech){
  case recursive:{
    METIS_PartGraphRecursive(&metis_nvtxs, metis_xadj, metis_adjncy, metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, options, &edgecut, metis_part);
    break;
  }
  case kway:{
    METIS_PartGraphKway(&metis_nvtxs, metis_xadj, metis_adjncy, metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, options, &edgecut, metis_part);
    //printf("edge cut %ld\n",edgecut);
    break;
  } 
  case vkway:{
    METIS_PartGraphVKway(&metis_nvtxs, metis_xadj, metis_adjncy, metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, options, &volume, metis_part);
    printf("volume %d\n",volume);
    break;
  }
  case mcrekursive:{
    METIS_mCPartGraphRecursive(&metis_nvtxs,&ncon, metis_xadj, metis_adjncy, metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, options, &edgecut, metis_part);
    break;
  }
  case mckway:{
    METIS_mCPartGraphKway(&metis_nvtxs,&ncon, metis_xadj, metis_adjncy, metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, ubec,options, &edgecut, metis_part);
    break;
  }
  case wpartrekursive:{
    METIS_WPartGraphRecursive(&metis_nvtxs, metis_xadj, metis_adjncy, metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, tpwgts,options, &edgecut, metis_part);
    break;
  }
  case wpartkway:{
    METIS_WPartGraphKway(&metis_nvtxs, metis_xadj, metis_adjncy,metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, tpwgts,options, &edgecut, metis_part);
    break;
  }
  case wpartvkway:{
    METIS_WPartGraphVKway(&metis_nvtxs, metis_xadj, metis_adjncy, metis_vwgt, metis_adjwgt, &wgtflag, &numflag, &nparts, tpwgts,options, &volume, metis_part);
    break;
  }
  }
  te = clock();
  printf("Time of metis partitoning %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  
  delete []metis_xadj;
  delete []metis_adjncy;
  delete []metis_vwgt;
  delete []metis_adjwgt;
    
  part = new int[nvtxs];
  for(i = 0; i < nvtxs; i++){
    part[i] = metis_part[i];
  }
  delete []metis_part;
  aggregate_postprocessing();
  
}


int meshdecomp::run()
{
  
  //long i,l;
  
  Gtm = new gtopology;
  Top = new siftop;
  Check_unused_prop = 1;
  
  if (Argc < 2){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use :%s input_file_name \n\n",Argv[0]);
    delete Gtm; delete Top;
    return(1);
  }
  
  err=read_input_data();

  err=read_topology();
  if(processing != aggregates){
    switch(partType){
    case METIS:{
      run_metis();
      break;
    }
    case CHACO:{
      run_chaco();
      break;
    }
    case JOSTLE:{
      run_jostle();
      break;
    }
    case PARTY:{
      run_party();
      break;
    }
    case SCOTCH:{
      run_scotch();
      break;
    }
    }
  }
  else{
    run_aggregates();
  }
  
  fclose(Log);

  return(0);
}








int main (int argc,char *argv[])
{
  if (argc < 2){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use :%s input_file_name \n\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  
  
  meshdecomp *decomp;
  double tws,twe;
  struct tm *local;
  time_t t;
  
  set_prgname(argv[0]);
  printf("--------------------------\n");
  printf("*** MESH DECOMPOSITION ***\n");
  printf("--------------------------\n");
  t = time(NULL);
  local = localtime(&t);
  printf("Program started  %s\n", asctime(local));
  // printf("with this input line parametrs\n",argv[0]);
  //   for(i = 1; i < argc; i++){
  //     printf(" %s",argv[i]);
  //   }
  //   printf("\n");
  
  tws = clock();
  decomp = new meshdecomp(argc,argv);
  decomp->run();
  decomp->~meshdecomp();
  twe = clock();
  printf("--------------------------\n");
  printf ("Whole elapsed time of mesh partitioning %le s\n",(twe-tws)/(double)CLOCKS_PER_SEC);
}
 
