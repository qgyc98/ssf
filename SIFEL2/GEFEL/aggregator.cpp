#include "aggregator.h"
#include <stdlib.h>
#include <limits.h>

aggregator::aggregator ()
{
  //  number of aggregates
  na=0;
  //  number of nodes
  nn=0;
  //  number of unknowns in the whole problem
  n=0;
  //  size of the coarse matrix
  cms=0;
  
  //  type of BOSS algorithm
  impl=0;
  
  //  type of rigid body modes
  trbm=0;
  //  number of rigid body modes
  nrbm=0;
  
  //  estimate of the spectral radius
  specrad=0.0;
  
  // degree of smoothing aggregator
  degree_k=0;
  degree_r=0;
  
  //  exact solver is used 
  exinex=1;
  
  
  lnntagr=NULL;
  lntagr=NULL;

  lnnagr=NULL;
  lnagr=NULL;

  //  list of numbers of unknowns on aggregates
  lnuaggr=NULL;
  //  list of unknown numbers on aggregates
  luaggr=NULL;
  
  //  matrices
  lmdm=NULL;
  lmcr=NULL;
  lmsky=NULL;
  sdirect=NULL;

  //  solver of linear equations (on each aggregate)
  ssle = new slesolv();

}

aggregator::~aggregator ()
{
  delete [] lnuaggr;
  delete [] luaggr;
  
  delete [] lmdm;
  delete [] lmcr;
  delete [] lmsky;
  delete [] sdirect;
}

/**
   function reads data about aggregation
   
   @param in - input stream
   
   31.5.2008, JK
*/
void aggregator::read (gtopology *gt,XFILE *in,long mespr)
{
  //  type of BOSS implementation
  //  impl=1 - own implementation of the decomposition
  //  impl=2 - decomposition provided by METIS
  xfscanf (in,"%ld",&impl);
  if (impl!=1 && impl!=2)
    print_err("wrong type of BOSS implemenation is required",__FILE__,__LINE__,__func__);
  
  //  type of solver - exact or inexact
  //  exinex=1 - exact solver in coarse problem
  //  exinex=2 - inexact solver in coarse problem
  xfscanf (in,"%ld",&exinex);
  if (exinex!=1 && exinex!=2)
    print_err("wrong type of solver (exact/inexact) is required",__FILE__,__LINE__,__func__);
  
  
  //  degree of recursion
  xfscanf (in,"%ld",&degree_k);
  if (mespr==1)
    fprintf (stdout,"\n degree of recursion in BOSS is %ld",degree_k);
  
  //  type of rigid body modes
  //  trbm=1 - one rigid body mode (heat transfer)
  //  trbm=2 - three rigid body modes (plane stress)
  xfscanf (in,"%ld",&trbm);
  
  
  //  data about solver of system of linear equations
  //gtopology *gt;
  ssle->read (gt,in,mespr);
  
  /*
  //  type of solver of system of linear equations
  xfscanf (in,"%k%m","typelinsol",&linsolvertype_kwdset,(int*)&tlinsol);
  
  switch (tlinsol){
  case gauss_elim:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by gaussian elimination method");
    break;
  }
  case ldl:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by LDL method");
    break;
  }
  case lu:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by LU method");
    break;
  }
  case ll:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by LL method (Cholesky)");
    break;
  }
  case spdirldl:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by sparse direct method (LDL)");
    break;
  }
  case spdirlu:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by sparse direct method (LU)");
    break;
  }
  case spdirll:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by sparse direct method (LL)");
    break;
  }
    
  default:{
    print_err("unknown type of solver of linear equations is required",__FILE__,__LINE__,__func__);
  }
  }
  */
}

/**
   function prints data about aggregation
   
   @param out - output stream
   
   31.5.2008, JK
*/
void aggregator::print (FILE *out)
{
  //  type of BOSS implementation
  fprintf (out,"\n%ld",impl);
  
  //  type of solver - exact or inexact
  fprintf (out,"\n%ld",exinex);

  //  degree of recursion
  fprintf (out,"\n%ld",degree_k);

  //  type of rigid body modes
  fprintf (out,"\n%ld",trbm);

  //  type of solver of system of linear equations
  //fprintf (out,"\n%d",tlinsol);

  fprintf (out,"\n");
  ssle->print (out);
  
}




/**
   function defines required number of aggregates
   
   @param i - required number of aggregates
   
   JK, 4.2.2007
*/
void aggregator::define_na (long i)
{
  na=i;
}

/**
   function defines graph of the %matrix
   
   @param nnod - number of nodes
   @param nneigh - numbers of adjacent nodes
   nneigh[i]=j - there are j adjacent nodes to the i-th node
   @param listneigh - list of node numbers of adjacent nodes
   listneigh[i][j]=k - the j-th adjacent node to the i-th node has number k
   
   JK, 4.2.2007
*/
void aggregator::define_graph (long nnod, long *nneigh, long **listneigh)
{
  //  number of nodes in the finite lement mesh
  nn = nnod;
  
  nadjnodnod = nneigh;
  
  adjnodnod = listneigh;

}

/**
   functions creates disjoint covering of node graph
   set both list of sizes and list of lists
   
   PM
*/
void aggregator::prepare_tlnagr(void)
{
  //printf("zacinam  vytvaret pokryti\n");
// odhadnu pocet prvku v agregatu a trochu to zmensim
  long i,j,k;
  long local_n = 10*nn/na/11;
  long cnt;
  long *lnodes;
  long *vrstva;
  long vrstva_start, vrstva_konec, node, neig;
  
  lnodes=NULL;
  vrstva=NULL;

  // vytvorim predbezne pokryti
  lnodes = new long [nn];
  memset (lnodes,0,sizeof(*lnodes)*nn);
  if ( lnodes == NULL) 
    {
      printf(" nedostatek pameti\n");
      exit(-1);
    }
  vrstva = new long [nn+1];
  memset (vrstva,0,sizeof(*vrstva)*(nn+1));
  //vrstva = (long *)calloc(nn+1, sizeof(long));
  if ( vrstva == NULL) 
    {
      printf(" nedostatek pameti\n");
      exit(-1);
    }
  //printf("alokovano\n");
  for (i=1;i<nn;i++) lnodes[i]=-1;
  for (j=0;j<na;j++)
    {
      //printf("zacinam cyklovat j=%ld\n",j);
      // najdi prvni volnej
      for( i=1;i<nn;i++)
	if (lnodes[i]==-1) break;
      cnt = 1; // pocitam pridelene prvky
      lnodes[i] = j;
      vrstva[1] = i;
      vrstva_start = 1;  //aktualne pridana vrstva
      vrstva_konec = 1;
      //printf(" dil 2\n");
      for(;cnt<local_n;)
	{
	  // projdi vrstvu a pridej neoznacene sousedy
	  //printf("dil 3\n");
	  for(i=vrstva_start; i<=vrstva_konec;i++)
	    {
	      node = vrstva[i];
	      //printf("node=%ld\n",node);
	      for(k=0;k<nadjnodnod[node];k++)
		{
		  neig = adjnodnod[node][k];
		  if (lnodes[neig]==-1)
		    {
		      cnt++;
		      lnodes[neig]=j;
		      vrstva[cnt]=neig;
		      if (cnt==local_n) break;
		    }
		}
	      if(cnt==local_n) break;
	    }
	  vrstva_start = vrstva_konec+1;
	  vrstva_konec = cnt;
	  //printf("cnt=%ld local_n=%ld\n",cnt,local_n);
	}
      //printf("cyklus ukoncen\n");
    }
  // doplnim zbytky
  for(;;)
    {
      //printf("pridavam\n");
      cnt=1;
      for(node=0;node<nn;node++)
	{
	  for(i=0;i<nadjnodnod[node];i++)
	    {
	      neig = adjnodnod[node][i];
	      if (lnodes[neig]==-1)
		{
		  lnodes[neig]=lnodes[i];
		  cnt = 0;
		}
	    }
	}
      if(cnt==1) break;
    }
  //for (i=0;i<n;i++) printf("%ld\n",lnodes[i]);
  // nasypu vysledky - lnodes obsahuje vse potrebne
  //lnntagr = (long*)calloc(na, sizeof(long));
  lnntagr = new long [na];
  memset (lnntagr,0,sizeof(*lnntagr)*na);

  //lntagr  = (long **)calloc(na, sizeof(long*));
  lntagr  = new long * [na];
  lntagr[0] = new long [nn];
  for (i=0;i<nn;i++){
    lntagr[0][i]=0;
  }
  
  //lntagr[0] = (long *)calloc(nn,sizeof(long));

  // otestovat pamet
  if ((lnntagr==NULL) ||(lntagr==NULL) ||(lntagr[0]==NULL)) 
    {
      printf("nedostatek pameti\n");
      exit(-1);
    }
  // OK
  for(i=0;i<na;i++) lnntagr[i]=0;
  for(i=0;i<nn;i++)  lnntagr[lnodes[i]]++;
  //for(i=0;i<na;i++) printf("%ld\n",lnntagr[i]);
  // naplnim seznamy sousedu
  cnt=0;
  for(i=0;i<na;i++)
    {
      for(j=0;j<nn;j++)
	if(lnodes[j]==i) 
	  {
	    lntagr[0][cnt]=j;
	    cnt++;
	  }
    }
  // presypat zacatky vektoru
  for(i=1;i<na;i++)
    lntagr[i]=lntagr[i-1]+lnntagr[i-1];
  
  delete [] lnodes;
  delete [] vrstva;
  //free(lnodes);
  //free(vrstva);
  //printf("pokryti vytvoreno\n");
}

/**
   function copies arrays lnntagr and lntagr from the
   object stop (class seqtop) allocated in class gtopology

   it is an alternative to the function void prepare_tlnagr(void);
   
   @param gt - global topology
   
   JK, 31.5.2008
*/
void aggregator::metis_aggr (gtopology *gt)
{
  long i,j;
  
  //  list of numbers of nodes in aggregates
  lnntagr = new long [na];
  for (i=0;i<na;i++){
    lnntagr[i]=gt->stop->nnsd[i];
  }
  
  //  list of node numbers in aggregates
  lntagr = new long* [na];
  for (i=0;i<na;i++){
    lntagr[i] = new long [lnntagr[i]];
  }
  for (i=0;i<na;i++){
    for (j=0;j<lnntagr[i];j++){
      lntagr[i][j]=gt->stop->ltg[i][j];
    }
  }
  
}

/**
   pripravi seznamy smoothed prolongatoru
   
   pred zavolanim teto funkce se musi zavolat funkce define_degree(long dg);
   
   PM
*/
void aggregator::prepare_lnagr ()
{
  long i,j,k,l,vrstva_zacatek,vrstva_konec,cnt;
  long node, neig;
  long *vrstva,*lnodes;
  vrstva=NULL;
  lnodes=NULL;
  // zname jeho stupen, tedy jen pridame vrstvy
  //vrstva = (long*)calloc(nn,sizeof(long));
  vrstva = new long [nn];
  memset (vrstva,0,sizeof(*vrstva)*nn);
  //lnodes = (long*)calloc(nn,sizeof(long));
  lnodes = new long [nn];
  memset (lnodes,0,sizeof(*lnodes)*nn);
  //lnnagr = (long*)calloc(na, sizeof(long));
  lnnagr = new long [na];
  memset (lnnagr,0,sizeof(*lnnagr)*na);
  //lnagr = (long**)calloc(na, sizeof(long*));
  lnagr = new long* [na];
  
  for(i=0;i<na;i++)
  {
    // napred nakopirujeme odpovidajici tentativni prolongator
    for(j=0; j<nn; j++) lnodes[j]=-1;
    for(j=0;j<lnntagr[i];j++) 
    {
      vrstva[j]=lntagr[i][j];
      lnodes[vrstva[j]] = i;
    }
    
    vrstva_zacatek = 0;
    vrstva_konec = lnntagr[i];
    for(k=0;k<degree_r;k++)
    {
      cnt = vrstva_konec;
      for(j=vrstva_zacatek;j<vrstva_konec;j++)
      {
        //projdi sousedy a oznac
        node=vrstva[j];
        for(l=0;l<nadjnodnod[node];l++)
        {
          neig = adjnodnod[node][l];
          if (lnodes[neig]==-1)
          {
            lnodes[neig] = i;
            vrstva[cnt] = neig;
            cnt++;
          }
        }
      } 
      vrstva_zacatek = vrstva_konec;
      vrstva_konec = cnt;
    }
    // ted mam ve vrstva vsechno co je treba
    lnnagr[i] = vrstva_konec;
    //lnagr[i] = (long*)calloc(nn,sizeof(long));
    lnagr[i] = new long [nn];
    for(l=0;l<vrstva_konec;l++)
      lnagr[i][l]=vrstva[l];
  }
  
  delete [] vrstva;
  delete [] lnodes;
}



/** 
    function defines depth of recursion for smoothed
    prolongator
    dg = 0 .... use polynomial of degree 0
    
    PM
*/
void aggregator::define_degree()
{
  degree_r = 1;
  for(long i=0;i<degree_k;i++) degree_r *= 3;
  degree_r = (degree_r-1)/2;
}


void aggregator::assemble_smoothed_prol2 (long /*agrid*/,long /*cid*/,long */*i*/,compvect *cvi,compvect *cvo,gmatrix *mtx)
{
  compvect *cva;

  cva = new compvect;

  //  prototyp nasobeni
  mtx->cr->crxcv_cv (cvi,cvo);
  
  




  
  //  kontrola na konec
  if (cvi->nz != cvo->nz){
    fprintf (stderr,"\n\n ruzne pocty nenulovych prvku ve vektorech cvi a cvo ve funkci assemble_smoothed_prol2");
    fprintf (stderr,"\n (file %s, line %d)\n",__FILE__,__LINE__);
  }
  
  delete [] cva;
  cva=NULL;
}


/** 
    vytvori zhlazeny prolongator
    input:
    agrid . cislo pozadovaneho agregatu
    cid ... cislo pozadovaneho sloupce v agregatu ?
    i   ... pole indexu ?
    a   ... pole hodnot (v pripade rigid body motions v mechanice)
    mtx ... matice problemu     
    PM 
*/
void aggregator::assemble_smoothed_prol (gtopology *gt,gmatrix *mtx)
{ 
  long i,j;
  double *v; // tohle je vektor se shlazenym prolongatorem
  // ted ho spravne vytvorit a naplnit
  v = new double[n];
  
  p = new double [na*nrbm*n];
  
  switch (trbm){
  case 1:{
    //  heat transfer, 1D mechanical problem

    for (i=0;i<na;i++){
      gener_rbm_1 (gt,i,v);
      
      // konstrukce probehne rekurzivne
      mulS(degree_k,v,mtx);
      
      copyv (v,p+i*n,n);
    }
    break;
  }
  case 2:{
    //  plane stress
    for (i=0;i<na;i++){
      for (j=0;j<nrbm;j++){

	gener_rbm_2 (gt,i,j,v);
	
	// konstrukce probehne rekurzivne
	mulS(degree_k,v,mtx);
	
	copyv (v,p+(i*nrbm+j)*n,n);
      }
    }
    break;
  }
  default:{
    print_err("unknown type of rigid body modes is required",__FILE__,__LINE__,__func__);
  }
  }
  
  delete [] v;
}

// lokalni metody - casem asi private 
// !!! uzitecne jen pro pripravu shlazeneho prolongatoru
/**
   provedena aplikaci s(st) na pozadovany vekor
   st    .... stupen S, t.j. dolni index
   x     .... vektor
   mtx   .... matice
**/
void aggregator::mulS(long st,double *x, gmatrix *mtx)
{
  double locro = specrad;
  double *wrk1;
  double *wrk2;
  
  wrk1 = new double[n];
  wrk2 = new double[n];
  
  for(long j=1;j<st;j++)
    {
      copyv(x,wrk1,n);     // wrk1 := x
      copyv(x,wrk2,n);     // wrk2 := x
      mulA(j-1,wrk2,mtx);  // wrk2 := A_j-1 * wrk2
      for (long k=0; k<n;k++) x[k] = wrk1[k]-(4*wrk2[k])/(3*locro);
      // x:= (I-(4/(3*locro))*A_j-1 * x
      locro /= 9;
    }
  delete [] wrk1;
  delete [] wrk2;
}


/** 
    st    .... stupen S, t.j. dolni index
    x     .... vektor
    mtx   .... matice
**/
void aggregator::mulA(long st, double *x, gmatrix *mtx)
{
  double *wrk1;
  double *wrk2;
  double locro = specrad;
  
  wrk1 = new double[n];
  wrk2 = new double[n];
  
  for(long k=1; k<st;k++) locro /= 9;
  copyv(x,wrk1,n);  // wrk1 := x
  if ( st == 0 )
    mtx->gmxv(wrk1,x);  // x := A*x
  else
    { 
      mulA(st-1,x,mtx); // x := A_st-1 * x
      copyv(x,wrk1,n);     // wrk1 := x
      copyv(x,wrk2,n);     // wrk2 := x
      mulA(st-1,wrk2,mtx);  // wrk2 := A_st-1 * wrk2
      for (long k=0; k<n;k++) x[k] = wrk1[k]-(4*wrk2[k])/(3*locro);
      // x := (I-4/(3*ro_st-1) * A_st-1 )A_st-1 *x 
      copyv(x,wrk1,n);     // wrk1 := x
      copyv(x,wrk2,n);     // wrk2 := x
      mulA(st-1,wrk2,mtx);  // wrk2 := A_st-1 * wrk2
      for (long k=0; k<n;k++) x[k] = wrk1[k]-(4*wrk2[k])/(3*locro);
      // x := (I-4/(3*ro_st-1) * A_st-1 )(I-4/(3*ro_st-1) * A_st-1 )A_j_1 *x 
    }
  delete [] wrk1;
  delete [] wrk2;
}


/**
   function computes kernels for heat transfer, 1D mechancs, etc.
   
   @param gt - pointer to the general topology
   @param agid - aggregate id
   @param v - array containing kernel
   
   JK, 27.8.2008
*/
void aggregator::gener_rbm_1 (gtopology *gt,long agid,double *v)
{
  long i,cn,nid;
  
  for (i=0;i<lnnagr[agid];i++){
    nid=lnagr[agid][i];
    cn=gt->give_dof (nid,0);
    if (cn>0)
      v[cn-1]=1.0;
  }
}

/**
   function computes kernels for plane stress
   
   @param gt - pointer to the general topology
   @param agid - aggregate id
   @param cnid - code number id
   @param v - array containing kernel
   
   JK, 27.8.2008
*/
void aggregator::gener_rbm_2 (gtopology *gt,long agid,long cnid,double *v)
{
  long i,cn,nid;
  double xp,yp,xk,yk;
  
  if (cnid==2){
    nid=lnagr[agid][0];
    xp=gt->gnodes[nid].x;
    yp=gt->gnodes[nid].y;
  }

  for (i=0;i<lnnagr[agid];i++){
    nid=lnagr[agid][i];
    if (cnid==0){
      cn=gt->give_dof (nid,cnid);
      if (cn>0)
	v[cn-1]=1.0;
    }
    if (cnid==1){
      cn=gt->give_dof (nid,cnid);
      if (cn>0)
	v[cn-1]=1.0;
    }
    if (cnid==2){
      xk=gt->gnodes[nid].x;
      yk=gt->gnodes[nid].y;
      
      cn=gt->give_dof (nid,0);
      if (cn>0)
	v[cn-1]=yp-yk;
      
      cn=gt->give_dof (nid,1);
      if (cn>0)
	v[cn-1]=xk-xp;
    }
  }
}

/**
   function computes coarse %matrix of the problem
   
   at this time, first and unefficient implementation is used
   the %matrix is stored in the dense %matrix storage scheme
   
   @param gm - general %matrix
   
   JK, 27.8.2008
*/
void aggregator::coarse_matrix (gmatrix *gm)
{
  long i,j;
  double s,*am;
  
  //  coarse matrix
  cm = new densemat [1];
  cm->alloc (cms);

  //  auxiliary matrix
  am = new double [n*cms];
  
  for (i=0;i<cms;i++){
    gm->gmxv (p+i*n,am+i*n);
  }
  
  for (i=0;i<cms;i++){
    for (j=0;j<cms;j++){
      s=ss (p+i*n,am+j*n,n);
      cm->add_entry (s,i,j);
    }
  }
  
  delete [] am;
}


/**
   function assembles list of unknowns belonging to the aggregates with overlap
   
   function assembles array luaggr which contains unknow numbers in usual (not C) notation
   it means, it starts from 1 and not from 0

   @param gt - pointer to general topology
   
   JK, 19.2.2007
*/
void aggregator::assemble_aggr_unknowns (gtopology *gt)
{
  long i,j,k,l,ann,nu,ndofn;
  
  if (lnuaggr!=NULL){
    delete [] lnuaggr;
  }
  lnuaggr = new long [na];
  
  maxnu=0;
  for (i=0;i<na;i++){
    nu=0;
    for (j=0;j<lnnagr[i];j++){
      //  actual node number
      ann=lnagr[i][j];
      //  number of unknowns on node
      ndofn = gt->give_ndofn (ann);
      for (k=0;k<ndofn;k++){
	l=gt->give_dof (ann,k);
	if (l>0)
	  nu++;
      }
      
    }
    lnuaggr[i]=nu;
    if (nu>maxnu)
      maxnu=nu;
  }
  
  if (luaggr!=NULL){
    for (i=0;i<na;i++){
      delete [] luaggr[i];
    }
    delete [] luaggr;
  }
  luaggr = new long* [na];
  for (i=0;i<na;i++){
    luaggr[i] = new long [lnuaggr[i]];
  }
  
  
  for (i=0;i<na;i++){
    lnuaggr[i]=0;
    for (j=0;j<lnnagr[i];j++){
      //  actual node number
      ann=lnagr[i][j];
      //  number of unknowns on node
      ndofn = gt->give_ndofn (ann);
      for (k=0;k<ndofn;k++){
	l=gt->give_dof (ann,k);
	if (l>0){
	  luaggr[i][lnuaggr[i]]=l;
	  lnuaggr[i]++;
	}
      }
      
    }
  }

}

/**
   function assembles local matrices of aggregates

   @param bsize - size of blocks for sparse direct solvers only
   @param gm - %matrix of the system
   
   JK, 20.2.2007
*/
void aggregator::local_matrices (long /*bsize*/,gmatrix *gm)
{
  long i,nu,pam,min,max;
  
  /*
  if (lmcr!=NULL){
    delete [] lmcr;
  }
  */
  lmcr = new comprow [na];

  /*
  if (lmdm!=NULL){
    delete [] lmdm;
  }
  lmdm = new densemat [na];
  */
  
  /*
  if (sdirect!=NULL){
    delete [] sdirect;
  }
  sdirect = new DSSolver [na];
  */
  
  pam=0;
  min=LONG_MAX;
  max=0;
  for (i=0;i<na;i++){
    gm->cr->select_submatrix (luaggr[i],lnuaggr[i],&lmcr[i]);
    //gm->cr->select_submatrix (luaggr[i],lnuaggr[i],&lmdm[i]);
    
    //  number of unknowns in the i-th aggregate
    nu=lmcr[i].n;
    if (lmcr[i].adr[nu]<min)
      min=lmcr[i].adr[nu];
    if (lmcr[i].adr[nu]>max)
      max=lmcr[i].adr[nu];
    
    pam+=lmcr[i].adr[nu];
    fprintf (stdout,"\n number of nonzero matrix entries in compressed rows storage %ld",lmcr[i].adr[nu]);
  }
  
  fprintf (stdout,"\n pocet prvku v maticich comp rows  %10ld",pam);
  fprintf (stdout,"\n minimalni pocet prvku             %10ld",min);
  fprintf (stdout,"\n maximalni pocet prvku             %10ld",max);
  
  
  if (exinex==1){
    //  exact solver is used
    lmsky = new skyline [na];
    for (i=0;i<na;i++){
      lmsky[i].assemble_from_cr (lmcr[i].n,lmcr[i].adr,lmcr[i].ci,lmcr[i].a);
    }
    
    pam=0;
    min=LONG_MAX;
    max=0;
    for (i=0;i<na;i++){
      pam+=lmsky[i].negm;
      if (min>lmsky[i].negm)
	min=lmsky[i].negm;
      if (max<lmsky[i].negm)
	max=lmsky[i].negm;
    }
    fprintf (stdout,"\n pocet prvku v maticich skyline  %10ld",pam);
    fprintf (stdout,"\n minimalni pocet prvku           %10ld",min);
    fprintf (stdout,"\n maximalni pocet prvku           %10ld",max);
  }
  
  
  // ***********************
  //  sparse direct solver
  // ***********************
  /*
  for (i=0;i<na;i++){
    
    switch (tlinsol){
    case spdirldl:{
      sdirect[i].Initialize(1,eDSSFactorizationLDLT);
      break;
    }
    case spdirlu:{
      sdirect[i].Initialize(1,eDSSFactorizationLU);
      break;
    }
    case spdirll:{
      sdirect[i].Initialize(1,eDSSFactorizationLLT);
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown type of sparse solver is required");
      fprintf (stderr,"\n in function aggregator::local_matrices (file %s, line %d)\n",__FILE__,__LINE__);
    }
    }
    
    sdirect[i].LoadMatrix ((unsigned long)lmcr[i].n,bsize,lmcr[i].a,(unsigned long*)lmcr[i].ci,(unsigned long*)lmcr[i].adr);
    //auxdatsparsesolver (top,out);
    
    //if (top->nn != 0){
      //sdirect->LoadMCN (top->nn,bsize,auxsd);
    //}

    
    SparseMatrixF sm((unsigned long)lmcr[i].n,lmcr[i].a,(unsigned long*)lmcr[i].ci,(unsigned long*)lmcr[i].adr);
    

      FILE *out;
      out = fopen ("matice.sm","wb");
      sm.SaveMatrix(out);
      fclose (out);

    
  }
  */  
  
  
  
  
  min=LONG_MAX;
  max=0;
  for (i=0;i<na;i++){
    if (min>lnuaggr[i])
      min=lnuaggr[i];
    if (max<lnuaggr[i])
      max=lnuaggr[i];
  }
  fprintf (stdout,"\n minimalni pocet neznamych v agregatu   %10ld",min);
  fprintf (stdout,"\n maximalni pocet neznamych v agregatu   %10ld",max);
  
}

void aggregator::clean_memory ()
{
  long i;
  
  for (i=0;i<na;i++){
    delete [] lntagr[i];
  }
  delete [] lntagr;
  
  delete [] lnntagr;
  
  
  for (i=0;i<na;i++){
    delete [] lnagr[i];
  }
  delete [] lnagr;
  
  delete [] lnnagr;
}


/**
   function prepares all necessary variables, arrays, etc. for BOSS algorithm
   
   @param gt - general topology
   @param gm - %matrix of the system stored in the class gmatrix
   @param out - output stream
   
   JK, 25.2.2007
*/
void aggregator::prepare_boss (gtopology *gt,gmatrix *gm,FILE *out)
{
  long i,j,k,l,v,min,bsize;
  time_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11;
  double *lhs,*rhs;
  lhs=NULL;
  rhs=NULL;

  //  number of aggregates
  na = gt->ns;
  
  //  number of nodes in the finite element mesh
  nn = gt->nn;
  
  //  number of unknowns in the whole problem
  n = gm->n;
  
  //  number of rigid body modes
  switch (trbm){
  case 1:{
    //  heat transfer, 1D mechanical problem
    nrbm=1;
    break;
  }
  case 2:{
    //  plane stress
    nrbm=3;
    break;
  }
  default:{
    print_err("unknown type of rigid body modes is required",__FILE__,__LINE__,__func__);
  }
  }
  
  //  size of the coarse matrix
  cms=na*nrbm;
  
  //  estimate of the spectral radius
  specrad = gm->estim_spect_radius ();
  
  t1 = time (NULL);
  
  fprintf (stdout,"\n searching for adjacent nodes");
  gt->adjacnodes (out);
  
  //  list of numbers of adjacent nodes
  nadjnodnod = gt->nadjnodnod;
  
  //  list of node numbers adjacent to nodes
  adjnodnod = gt->adjnodnod;
  
  t2 = time (NULL);
  
  // ****************************************
  //  definition of the graph of the matrix
  // ****************************************
  
  if (impl==1){
    fprintf (stdout,"\n prepare_tlnagr");
    
    //  list of node numbers belonging to particular aggregates without overlap
    //  arrays lnntagr and lntagr are assembled
    prepare_tlnagr ();
  }
  if (impl==2){
    //  aggregation is obtained from METIS
    metis_aggr (gt);
  }
  
  t3 = time (NULL);

  //  function sets up required degree of recursion
  define_degree ();
  
  fprintf (stdout,"\n prepare_lnagr");
  
  //  list of node numbers belonging to particular aggregates with overlap
  //  arrays lnnagr and lnagr are assembled
  prepare_lnagr ();
  
  t4 = time (NULL);
  
  fprintf (stdout,"\n trideni");
  
  //
  //  zde je nutne setridit pole lnagr a lntagr
  //
  for (i=0;i<na;i++){
    for (j=0;j<lnnagr[i];j++){
      min=LONG_MAX;
      for (k=j;k<lnnagr[i];k++){
	if (min>lnagr[i][k]){
	  min=lnagr[i][k];
	  l=k;
	}
      }
      v=lnagr[i][j];
      lnagr[i][j]=lnagr[i][l];
      lnagr[i][l]=v;
    }
  }
  

  /*
  delete [] lnnagr;
  for (i=0;i<na;i++){
    delete [] lnagr[i];
  }
  delete [] lnagr;
  

  lnnagr = new long [2];

  lnnagr[0]=3;
  lnnagr[1]=3;

  lnagr = new long* [2];
  for (i=0;i<2;i++){
    lnagr[i]=new long [lnnagr[i]];
  }

  
  lnagr[0][0]=0;
  lnagr[0][1]=1;
  lnagr[0][2]=2;

  lnagr[1][0]=2;
  lnagr[1][1]=3;
  lnagr[1][2]=4;
  */



  /*
  delete [] lnnagr;
  for (i=0;i<na;i++){
    delete [] lnagr[i];
  }
  delete [] lnagr;
  
  lnnagr = new long [3];

  lnnagr[0]=9;
  lnnagr[1]=9;
  lnnagr[2]=9;

  lnagr = new long* [3];
  for (i=0;i<3;i++){
    lnagr[i]=new long [lnnagr[i]];
  }

  
  lnagr[0][0]=0;
  lnagr[0][1]=1;
  lnagr[0][2]=2;
  lnagr[0][3]=3;
  lnagr[0][4]=4;
  lnagr[0][5]=5;
  lnagr[0][6]=6;
  lnagr[0][7]=7;
  lnagr[0][8]=8;

  lnagr[1][0]=7;
  lnagr[1][1]=8;
  lnagr[1][2]=9;
  lnagr[1][3]=10;
  lnagr[1][4]=11;
  lnagr[1][5]=12;
  lnagr[1][6]=13;
  lnagr[1][7]=14;
  lnagr[1][8]=15;

  lnagr[2][0]=11;
  lnagr[2][1]=12;
  lnagr[2][2]=13;
  lnagr[2][3]=14;
  lnagr[2][4]=15;
  lnagr[2][5]=16;
  lnagr[2][6]=17;
  lnagr[2][7]=18;
  lnagr[2][8]=19;
  */

  /*
  //  pomocny tisk pro ladeni
  fprintf (out,"\n\n pocty uzlu v jednotlivych agregatech\n");
  for (i=0;i<na;i++){
    fprintf (out,"agregat %5ld  obsahuje %6ld uzlu\n",i,lnnagr[i]);
  }
  fprintf (out,"\n\n seznamy uzlu v jednotlivych agregatech\n");
  for (i=0;i<na;i++){
    fprintf (out,"\n agregat %ld\n",i);
    for (j=0;j<lnnagr[i];j++){
      fprintf (out,"%ld\n",lnagr[i][j]);
    }
  }
  //  konec pomocneho tisku
  */
  
  
  

  //  zde bude volani teto funkce nebo jeji modifikace
  //void aggregator::assemble_smoothed_prol (long agrid,long cid,long *i,double *a, gmatrix *mtx)
 
  
  
  fprintf (stdout,"\n assemble_aggr_unknowns");

  t5 = time (NULL);

  //  list of unknowns belonging to particular aggregates
  //  arrays lnuaggr and luaggr are assembled
  assemble_aggr_unknowns (gt);
  

  t6 = time (NULL);
  
  //  size of blocks for sparse direct solvers
  bsize=gt->gnodes[0].ndofn;
  
  //  definition and assembling of matrices of aggregates
  local_matrices (bsize,gm);
  

  /*
  fprintf (out,"\n\n matice \n");
  gm->printmat (out);
  fprintf (out,"\n\n submatice 1\n");
  lmcr[0].printmat (out);
  lmsky[0].printmat (out);
  //lmdm[0].printmat (out);
  fprintf (out,"\n\n submatice 2\n");
  lmcr[1].printmat (out);
  lmsky[1].printmat (out);
  //lmdm[1].printmat (out);
  //fprintf (out,"\n\n submatice 3\n");
  //lmdm[2].printmat (out);
  */

  t7 = time (NULL);
  
  if (exinex==1){
    //  exact solver is used
    
    //  factorization of the local matrices
    for (i=0;i<na;i++){
      //lmdm[i].ll (lhs,rhs,1.0e-10,2);
      lmsky[i].ldl_sky (lhs,rhs,1.0e-10,2);
    }
  }
  
  t8 = time (NULL);
  
  //  zde se musi provest AA = P^T A P
  //  dale se musi matice AA rozlozit
  
  assemble_smoothed_prol (gt,gm);
  
  t9 = time (NULL);

  coarse_matrix (gm);

  t10 = time (NULL);
  
  double zero=1.0e-15;
  cm->ll(lhs,rhs,zero,2);

  t11 = time (NULL);

  clean_memory ();


  fprintf (stdout,"\n\n Casy v BOSS");
  fprintf (stdout,"\n sestaveni sousednich uzlu  %ld", long(t2-t1));
  fprintf (stdout,"\n sestaveni agregatu         %ld", long(t3-t2));
  fprintf (stdout,"\n sestaveni agregatu         %ld", long(t4-t3));
  fprintf (stdout,"\n trideni                    %ld", long(t5-t4));
  fprintf (stdout,"\n sestaveni neznamych        %ld", long(t6-t5));
  fprintf (stdout,"\n lokalni matice             %ld", long(t7-t6));
  fprintf (stdout,"\n rozklad lokalnich matic    %ld", long(t8-t7));
  fprintf (stdout,"\n zhlazeny prolongator       %ld", long(t9-t8));
  fprintf (stdout,"\n sestaveni coarse matice    %ld", long(t10-t9));
  fprintf (stdout,"\n rozklad coarse matice      %ld", long(t11-t10));
  fprintf (stdout,"\n");

}


/**
   function computes Black-box Overlapping Schwarz with Smoothed coarse space
   see Brezina, page 75 and 76
   
   @param gm - %matrix of the system
   @param u - input %vector
   @param v - output %vector
   
   JK, 20.2.2007
*/
void aggregator::boss (gmatrix *gm,double *u,double *v)
{
  long i;
  double *d,*dd,*pp,*z,*w,*ww;
  
  //  auxiliary arrays
  d = new double [n];
  dd = new double [maxnu];
  pp = new double [maxnu];
  z = new double [n];
  w = new double [cms];
  ww = new double [cms];

  // *********
  //  step 1
  // *********
  //copyv (u,z,n);
  fillv (0.0,z,n);
  
  // *********
  //  step 2
  // *********
  for (i=0;i<na;i++){
    //  A z^{i-1}
    gm->gmxv (z,v);
    //  f - A z^{i-1}
    subv (u,v,d,n);
    
    //  short vector cleaning
    fillv (0.0,dd,maxnu);
    //  reduction of long vectors to short vectors
    globloc (d,dd,luaggr[i],lnuaggr[i]);
    
    //  back-substitution
    //??????????????
    //lmdm[i].ll (pp,dd,1.0e-10,3);
    
    if (exinex==1){
      //  exact solver is used
      lmsky[i].ldl_sky (pp,dd,1.0e-10,3);
    }
    if (exinex==2){
      //  inexact solver if used
      lmcr[i].cg (pp,dd,ssle->ni,ssle->res,ssle->ani,ssle->ares,ssle->zero,0);
    }
    
    //  prolongation of short vectors to long vectors
    locglob (z,pp,luaggr[i],lnuaggr[i]);
  }
  
  // *********
  //  step 3
  // *********
  
  //  A z^{i-1}
  gm->gmxv (z,v);
  //  f - A z^{i-1}
  subv (u,v,d,n);
  
  /*
  for (i=0;i<cms;i++){
    w[i]=ss(p+i*n,d,n);
  }
  
  double zero=1.0e-10;
  cm->ll (ww,w,zero,3);
  
  for (i=0;i<n;i++){
    for (j=0;j<cms;j++){
      z[i]+=p[j*n+i]*ww[j];
    }
  }
  */

  // *********
  //  step 4
  // *********
  
  //  v Brezinove praci chybi???
  
  // *********
  //  step 5
  // *********
  for (i=na-1;i>-1;i--){
    //  A z^{i-1}
    gm->gmxv (z,v);
    //  f - A z^{i-1}
    subv (u,v,d,n);
    
    //  short vector cleaning
    fillv (0.0,dd,maxnu);
    //  reduction of long vectors to short vectors
    globloc (d,dd,luaggr[i],lnuaggr[i]);
    
    //  back-substitution
    //??????????????
    //lmdm[i].ll (pp,dd,1.0e-10,3);
    
    if (exinex==1){
      //  exact solver is used
      lmsky[i].ldl_sky (pp,dd,1.0e-10,3);
    }
    if (exinex==2){
      //  inexact solver is used
      lmcr[i].cg (pp,dd,ssle->ni,ssle->res,ssle->ani,ssle->ares,ssle->zero,0);
    }
    
    //  prolongation of short vectors to long vectors
    locglob (z,pp,luaggr[i],lnuaggr[i]);
  }
  

  // *********
  //  step 6
  // *********
  copyv (z,v,n);
  
  
  delete [] z;
  delete [] pp;
  delete [] dd;
  delete [] d;
  delete [] ww;
  delete [] w;
}
