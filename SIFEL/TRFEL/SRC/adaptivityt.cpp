#include "adaptivityt.h"

#include "globalt.h"
#include "globmatt.h"
#include "elemswitcht.h"


/// CONSTRUCTOR
adaptivityt :: adaptivityt ()
{
  dim = ord = ncomp = nn = ne = ntm = 0;
  
  tad        = Tp->adaptivityflag;
  printflags = 3;
  corr       = 1.0;
  adapt_accuracy = 1.0;
  //otherflags = -1;

  path = filename = suffix = NULL;
  interadaloop = false;
  step = 0;
  niwidth = 2;
  ni = NULL;
    
  refined_ders_in_nodes = elem_error_pct = refsizelrel = refsizelabs = NULL;
  
  r = rdr = NULL;
  tctrl = NULL;
  
  answer = 0;
}
/// DESTRUCTOR
adaptivityt :: ~adaptivityt ()
{
  delete [] ni;
  delete [] path;
  delete [] filename;
  delete [] suffix;
  delete [] r;
  delete [] rdr;
  delete    tctrl;
}

/// setup step and ni
void adaptivityt :: set_step (long s)
{
  step = s;
  if (ni == NULL)  ni = new char[1+niwidth+1];
  if (step >= 0)  sprintf (ni, ".%0*ld", niwidth, step);
  else            ni[0] = '\0';
}

/// read values from input file
void adaptivityt :: readinit (XFILE *in)
{
  xfscanf (in, "%k%lf %k%lf %k%lf", "accuracy", &adapt_accuracy, "enlargement", &enlarg, "reduction", &reduct);
  
  if (Tp->adaptivityflag<1 || Tp->adaptivityflag>2)
    { print_err("unknown type of adaptivity is required",__FILE__,__LINE__,__func__); exit (0); }
  
  filename_decomposition (in->fname, path, filename, suffix);
}

/// print values to input file
void adaptivityt :: printinit (FILE *out)
{
  fprintf (out,"%ld %.3lf %lf %.2lf\n", tad, adapt_accuracy, enlarg, reduct);
}

///
void adaptivityt :: initialize (long s)
{
  /// INITIALIZATION of VARIABLES
  dim   = Tt->give_dimension (0);
  ord   = Tt->give_degree (0);
  ncomp = Tt->give_ncomp (0);
  nn    = Tt->nn;
  ne    = Tt->ne;
  ntm   = Tp->ntm;
  
  for (long i=0; i<ne; i++)
    Gtt->gelements[i].auxinf = 100*Tt->give_nne(i) + 10*ord + dim;
  
  answer = 0;
  
  step = s;
  
  prepare_ni ();
  
  /// CHECK
  this->check_consistency ();
}


/**
   This is main function of class adaptivity, it runs in three steps:
   1. preparing strain
   2. computing of refined derivative (strain or stress, determined by otherflag) by spr_smoothing or z2_smoothing
   3. computing of error of FEM solution
   (case 3 = computing of accurate error)
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
long adaptivityt :: run (long /*of*/, bool ial)
{
  long hod,min;
  double sec;
  
  if (Mesprt>0){
    fprintf (stdout,"\n\n *********************************");
    fprintf (stdout,  "\n ***   Adaptivity is started   ***");
    
    sec = clock();
  }
  
  interadaloop = ial;
  
  /// strainy jsou pocitany standardne nekde vejs, kduz se pocitaji nejaky forces
  //compute_ipfluxes ();  // toto je na nic, protoze na tr ani pocitani fluxu neni
  compute_ipgrads ();   // tato fce je zapnuta pro pomocne adapt pocitani po zmene site, kdy gradienty jeste nejsou spocitany
  //                      // ALE nekontroluje asi ze se grads budou pocitat znova !! chtelo by to asi neco jako  Mp->strainstate=1;
  
  ///
  if (dim == 0)
    this->initialize(0);
  
  /// ALLOCATION
  refined_ders_in_nodes = new vector[nn];
  for (long i=0; i<nn; i++)  allocv (ncomp*ntm, refined_ders_in_nodes[i]); // *ntm je tam jen pro moznost tisku
  
  elem_error_pct = new vector[ntm+1];
  refsizelrel    = new vector[ntm+1];
  refsizelabs    = new vector[ntm+1];
  
  /// MAIN
  for (int i=0; i<ntm; i++) {
    // toto je zde jen pro moznost tisku
    if (ntm==2 && i==1)
      for (long j=0; j<nn; j++)
	for (long k=0; k<ncomp; k++)
	  refined_ders_in_nodes[j].a[k+i*ncomp] = refined_ders_in_nodes[j].a[k];
    
    switch (tad){
    case 1:{
      spr (i);
      compute_error (i);
      break;
    }
    default: print_err("unknown type of error smoothing is required in function", __FILE__, __LINE__, __func__);  exit (1);
    }
    
    
  }
  
  /// POST PROC
  if (ntm > 1) {
    allocv (ne, elem_error_pct[ntm]);
    allocv (ne, refsizelrel[ntm]);
    allocv (ne, refsizelabs[ntm]);
    
    for (long i=0; i<ne; i++) {
      elem_error_pct[ntm][i] = elem_error_pct[0][i];
      refsizelrel   [ntm][i] = refsizelrel   [0][i];
      refsizelabs   [ntm][i] = refsizelabs   [0][i];
      
      if (ntm > 1) {
	if (elem_error_pct[ntm][i] < elem_error_pct[1][i])  elem_error_pct[ntm][i] = elem_error_pct[1][i];
	if (refsizelrel   [ntm][i] > refsizelrel   [1][i])  refsizelrel   [ntm][i] = refsizelrel   [1][i];
	if (refsizelabs   [ntm][i] > refsizelabs   [1][i])  refsizelabs   [ntm][i] = refsizelabs   [1][i];
      }
      
      if (ntm > 2) {
	print_err("number of transported matters greater then 2", __FILE__, __LINE__, __func__);  exit (1);
      }
    }
  }
  
  /// *** print of background mesh file ***
  if (printflags & 2  &&  answer){
    char name[255];
    sprintf (name,"%s%s.T3d.bgm%s", path, filename, ni);
    print_valel (Outt, Gtt, path, name, "chyba", refsizelabs[ntm==1 ? 0 : ntm].a, 'e');
  }
  
  this->print_addat_vtk ();
  
  
  delete [] refined_ders_in_nodes;
  delete [] elem_error_pct;
  delete [] refsizelrel;
  delete [] refsizelabs;
  
  if (Mesprt>0){
    sec = (clock() - sec) / (double)CLOCKS_PER_SEC;
    hod = (long)sec/3600;  sec -= hod*3600;
    min = (long)sec/60;    sec -= min*60;
    fprintf(stdout,"\n Consumed time by Ada   %2ld:%02ld:%05.2f",hod,min,sec);
    
    fprintf (stdout,"\n ***   Adaptivity is finished  ***");
    fprintf (stdout,"\n *********************************");
  }
  return (answer);
}


///
void adaptivityt :: prepare_ni (void)
{
  if (! interadaloop) {
    step = atol(suffix+1);
    
    if ( step || 
	 (suffix[1] == '0' &&                                         suffix[2] == '\0') ||
	 (suffix[1] == '0' && suffix[2] == '0' &&                     suffix[3] == '\0') ||
	 (suffix[1] == '0' && suffix[2] == '0' && suffix[3] == '0' && suffix[4] == '\0') )
      ;
    else
      step = -1;
  }
  
  this->set_step(this->step);
  
  /// termitovo speciality
  if (strncmp(filename, "this.sifel", 10) == 0)
    filename[4] = '\0';
}

/**
   Function checks consistency of the problem.
   
   created  3.4.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void adaptivityt :: check_consistency (void) const
{
  int te = Tt->give_elem_type (0);
  for (long i=0; i<nn; i++) {
    if (ntm != Tt->give_ndofn    (i)) { print_err("mismatch in ndofn, node %ld differs", __FILE__, __LINE__, __func__, i+1);  exit (1); }
  }
  for (long i=0; i<ne; i++) {
    if (te  != Tt->give_elem_type(i)) { print_err("Uniform type of elements is required, element %ld differs", __FILE__, __LINE__, __func__, i+1);  exit (1); }
    if (ord != Tt->give_degree   (i)) { print_err("mismatch in order", __FILE__, __LINE__, __func__);  exit (1); }
    if (dim != Tt->give_dimension(i)) { print_err("mismatch in dimension", __FILE__, __LINE__, __func__);  exit (1); }
  }
  
  switch (te){
  case trlint:
  case quadlint:
    if (Ltt && (Ltt->intordkm[0][0]!=1 || Ltt->ncomp!=2)) { print_err("Ltt has nonstandard intpoints", __FILE__, __LINE__, __func__);  exit (1); }
    if (Lqt && (Lqt->intordkm[0][0]!=2 || Lqt->ncomp!=2)) { print_err("Lqt has nonstandard intpoints", __FILE__, __LINE__, __func__);  exit (1); }
    break;
  default:  print_err("unknown element type is required", __FILE__, __LINE__, __func__);  exit (1);
  }
  
  //  udelat kontrolu
  // predpoklad - na elementu je konstantni matice tuhosti d
  // funguje pro domenu s jednim materialem
  //mam kontrolovat stejny material po domene?
  //jak? po integracnich bodech?
  //asi na ipp nebo kolik je materialu na domene
}

/**
   Function fills up array refined_ders_in_nodes by derivatives(strain or stresses - it depends on otherflags) in nodes.
   It fills up array 'spder' by strain or stresses in sampling points
   and calls spr_smoothing to compute refined derivatives in nodes.
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void adaptivityt :: spr (int mattid)
{
  // This loop prepares array 'spder' for using in function `least_square::solve`
  // It allocates and fills up array 'spder' by derivatives(strain or stress) in main Gauss's points.
  // There is no choice now, strain is prefered.
  matrix *spder = new matrix[ne];
  
  /// diky tomuto muzu vse brat z prvniho bloku
  if (Tp->savemode!=1) {  print_err("savemode is NOT TRUE", __FILE__, __LINE__, __func__);  exit (0); }
  
  int nip;
  long i, j, ipp;
  vector vect;
  vect.n = ncomp;
  for (i=0; i<ne; i++) {
    // nip = Mt->give_nip(i, 0, 0); // prvni blok
    // tato fce zatim neni, nahore se nip kontroluje natvrdo
    switch (Tt->elements[i].te) {
    case quadlint:  nip = 4;  break;
    case trlint:    nip = 1;  break;
    default: print_err("unknown element type", __FILE__, __LINE__, __func__);  exit (0);
    }
    
    allocm(nip, ncomp, spder[i]);
    ipp = Tt->elements[i].ipp[0][0];  // save mode supposed
    
    for (j=0; j<nip; j++) {
      vect.a = spder[i].a + j*ncomp;
      Tm->givegrad (mattid, ipp, vect);
      ipp++;
    }
  }
  vect.a = NULL;
  
  // least_square spr(dim,ncomp,Gtm,0);
  patch_averaging spr (Gtt, dim, ncomp, 0);
  
  spr.solve (Outt, spder, refined_ders_in_nodes);
  
  delete [] spder;
}

/**
   Function computes error of FEM solution from difference between
   apriximate derivative fied(from FEM solution) and much more exact derivative fied(from array refined_ders_in_nodes).
   Considering 'printflags' printing of any files and information is performed.
   
   created  3.6.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void adaptivityt :: compute_error (int mattid)
{
  long i;
  double e2     ,u2     ,volume   ;
  vector ei2(ne),ui2(ne),sizel(ne);
  
  e2 = u2 = volume = 0.0;
  
  for (i=0; i<ne; i++) {
    
    switch (Tt->give_elem_type (i)) {
    case trlint:    volume += Ltt->compute_error (i, refined_ders_in_nodes, mattid, ei2[i], ui2[i], sizel[i]);  break;
    case quadlint:  volume += Lqt->compute_error (i, refined_ders_in_nodes, mattid, ei2[i], ui2[i], sizel[i]);  break;
    default:  print_err("wrong mesh - element %ld", __FILE__, __LINE__, __func__, i+1);  exit (1);
    }
    
    e2 += ei2[i];
    u2 += ui2[i];
  }
  u2 += e2;
  
  //
  compute_refsizel_lin (mattid, sizel.a, ei2.a, e2, u2);
  
  /// percentual global error
  double er = sqrt(e2/u2);
  if (er > adapt_accuracy || er < 0.02)  answer = 1;
  //if (er > adapt_accuracy)  answer = 1;
  
  if (Mesprt>0)
    fprintf (stdout, "\n Matter %d: global error = %.1lf%% (sqrt(e2/u2) == sqrt(%.3e / %.3e)) at volume = %f", mattid+1, 100*er, e2, u2, volume);
  
  ///
  allocv (ne, elem_error_pct[mattid]);
  for (i=0; i<ne; i++)
    elem_error_pct[mattid][i] = 100.0 * sqrt( ei2[i] / (ui2[i]+ei2[i]) );
  
}


/**
   Function computes new=required characteristic size of elements.
   
   @param sizel - array of old characteristic sizes of element
   @param
   
   created  3.6.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void adaptivityt :: compute_refsizel_lin (int mattid, const double *sizel, const double *ei2, double /*e2*/, double u2)
{
  long i;
  double em2,specacc;
  vector epsil;
  
  allocv (ne, epsil);
  
  specacc = this->adapt_accuracy;
  em2 = specacc*specacc * u2/ne;
  
  corr = 1.0;  //<0.7;1.0>
  
  allocv (ne, refsizelrel[mattid]);
  allocv (ne, refsizelabs[mattid]);
  
  for (i=0; i<ne; i++) {
    epsil[i] = sqrt(ei2[i] / em2);
    if (epsil[i] < 1.0 / enlarg) epsil[i] = 1.0 / enlarg;
    if (epsil[i] > 1.0 / reduct) epsil[i] = 1.0 / reduct;
    
    refsizelrel[mattid][i] = 1.0      / pow(epsil[i] , 1.0/(double)ord * corr);
    refsizelabs[mattid][i] = sizel[i] * refsizelrel[mattid][i];
  }
  
  // ***********************************************
  //char file[255];
  //sprintf (file,"%sepsil.dx%s",Mp->path,ni);
  //print_valel (Gtm,Mp,file,"epsil",epsil.a,'d');
  // ***********************************************
}


/// print additional data in vtk format
void adaptivityt :: print_addat_vtk () const
{
  long i,j;
  char name[255];
  FILE *out;
  
  //sprintf (name,"%s%s.addat.vtk%s", path, filename, suffix);
  sprintf (name,"%s%s.addat.vtk%s", path, filename, ni);
  out = fopen (name,"w");   if (out==NULL) fprintf (stderr,"\n .test file has not been opened.\n\n");
  
  fprintf (out,"#\n");
  // ***************************************************************************
  fprintf (out,"POINT_DATA %ld\n", nn);
  
  /// refined derivatives
  fprintf (out,"SCALARS ada_refined_derivatives float %ld\n", ncomp*ntm);
  fprintf (out,"LOOKUP_TABLE default\n");
  
  for (i=0; i<nn; i++) {
    for (j=0; j<ncomp*ntm; j++)
      fprintf (out, " %13.6f", refined_ders_in_nodes[i].a[j]);
    
    fprintf (out,"\n");
  }
  
  
  // ***************************************************************************
  fprintf (out,"CELL_DATA %ld\n", ne);
  
  /// percentual error at element
  fprintf (out,"SCALARS ada_error_pct float %ld\n", ntm==1 ? 1 : ntm + 1);
  fprintf (out,"LOOKUP_TABLE default\n");
  
  for (i=0; i<ne; i++) {
    ;             fprintf (out, "%13.6f", elem_error_pct[0]  [i]);
    if (ntm > 1)  fprintf (out, "%13.6f", elem_error_pct[1]  [i]);
    if (ntm !=1)  fprintf (out, "%13.6f", elem_error_pct[ntm][i]);
    
    fprintf (out,"\n");
  }
  
  /// refined size of element
  fprintf (out,"SCALARS ada_refsizel_relative float %ld\n", ntm==1 ? 1 : ntm + 1);
  fprintf (out,"LOOKUP_TABLE default\n");
  
  for (i=0; i<ne; i++) {
    ;             fprintf (out, "%13.6f", refsizelrel[0]  [i]);
    if (ntm > 1)  fprintf (out, "%13.6f", refsizelrel[1]  [i]);
    if (ntm !=1)  fprintf (out, "%13.6f", refsizelrel[ntm][i]);
    
    fprintf (out,"\n");
  }
  
  fclose (out);
}


// ********************************************************************************************************
// ********************************************   STATE DATA   ********************************************
// ********************************************************************************************************


///   !!!! toto narvat do gtopology pod dohledem JKTK   !!!!!!!!!!!!!!!

/**
   Function finds closest node to 'point'.
   The point is defined by coordinates x,y.
   
   @param dim - dimension of solved problem
   @param x,y - coordinates of the 'point'
   @retval    - number of closest node
   
   created  11.12.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
long closest_node (gtopology *gt, long dim, double x, double y, double z)
{
  long i, ans;
  double aux, dist=1e100;
  double coord[3];
  
  coord[0] = x;
  coord[1] = y;
  coord[2] = z;
  
  for (i=0; i<gt->nn; i++)
    if ( (aux=gt->gnodes[i].distance2(dim, coord)) < dist ) {
      dist = aux;
      ans = i;
    }
  
  return ans;
}

/**
   Function finds out, whether point (xx,yy) lays on the left hand side of the line
   defined by two nodes n1 and n2.
   
   created  3.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
bool isnodonlhsofline (gnode &n1, gnode &n2, double xx, double yy)
{
  if (n1.x != n2.x)
    if ( (n2.x>n1.x ? 1.0:-1.0)*((xx-n1.x)*(n1.y-n2.y)/(n1.x-n2.x)+(n1.y-yy)) <= 1.0e-13 )   return true;
    else                                                                                     return false;
  else
    if ( (n1.y-n2.y)*(xx-n1.x) > -1.0e-13 )  return true;
    else                                     return false;
}

/**
   Function finds out, whether point (xx,yy) lays on the left hand side of the curve
   defined by three nodes n1,n2 and n3.
   
   created  3.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
bool isnodonlhsof3pcurve (gnode &n1, gnode &n2, gnode &n3, double xx, double yy)
{
  double c;
  
  c = (n1.x*(n3.y-n2.y)+n2.x*(n1.y-n3.y)+n3.x*(n2.y-n1.y))/(fabs(n1.x-n3.x)+fabs(n1.y-n3.y));
  if (-1.0e-13<c && c<1.0e-13)
    return (isnodonlhsofline (n1,n3,xx,yy));
  else {
    // for testing
    // fprintf (stdout,"\n\n\n hrana je krivocara \n\n\n");
    // for testing
    double s,x[3],y[2];
    
    c = sqrt((n3.x-n1.x)*(n3.x-n1.x)+(n3.y-n1.y)*(n3.y-n1.y));
    s = (n3.y-n1.y)/c;
    c = (n3.x-n1.x)/c;
    
    x[0] =  c*(  xx-n1.x) + s*(  yy-n1.y);  y[0] = -s*(  xx-n1.x) + c*(  yy-n1.y);
    x[1] =  c*(n2.x-n1.x) + s*(n2.y-n1.y);  y[1] = -s*(n2.x-n1.x) + c*(n2.y-n1.y);
    x[2] =  c*(n3.x-n1.x) + s*(n3.y-n1.y);
    
    if ( y[1]*x[0]*x[0]/x[1]/(x[1]-x[2]) + y[1]*x[2]*x[0]/x[1]/(x[2]-x[1]) -y[0] <= 1.0e-13 )
      return true;
    
    return false;
  }
}


/**
   Function finds out, whether 'node' lays in 'element'.
   Node is defined by its coordinates - x,y.
   Element is defined by its node numbers - 'elnod' and dimdegnne - 'ndd'
   
   @param nod - array of all gnodes
   @param elnod - array of node numbers of element
   @param x,y - coordinates of node
   @param ndd - nnedimdeg of element
   
   created  3.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
long ispointinel (gnode *nod, long *elnod, double x, double y, long ndd)
{
  switch (ndd) {
  case 312:{
    if ( isnodonlhsofline (nod[elnod[0]],nod[elnod[1]],x,y) &&
	 isnodonlhsofline (nod[elnod[1]],nod[elnod[2]],x,y) &&
	 isnodonlhsofline (nod[elnod[2]],nod[elnod[0]],x,y)     )
      return true;
    break;
  }
  case 412:{
    if ( isnodonlhsofline (nod[elnod[0]],nod[elnod[1]],x,y) &&
	 isnodonlhsofline (nod[elnod[1]],nod[elnod[2]],x,y) &&
	 isnodonlhsofline (nod[elnod[2]],nod[elnod[3]],x,y) &&
	 isnodonlhsofline (nod[elnod[3]],nod[elnod[0]],x,y)     )
      return true;
    break;
  }
  case 622:{
    if ( isnodonlhsof3pcurve (nod[elnod[0]],nod[elnod[3]],nod[elnod[1]],x,y) &&
	 isnodonlhsof3pcurve (nod[elnod[1]],nod[elnod[4]],nod[elnod[2]],x,y) &&
	 isnodonlhsof3pcurve (nod[elnod[2]],nod[elnod[5]],nod[elnod[0]],x,y)     )
      return true;
    break;
  }
  case 822:{
    if ( isnodonlhsof3pcurve (nod[elnod[0]],nod[elnod[4]],nod[elnod[1]],x,y) &&
	 isnodonlhsof3pcurve (nod[elnod[1]],nod[elnod[5]],nod[elnod[2]],x,y) &&
	 isnodonlhsof3pcurve (nod[elnod[2]],nod[elnod[6]],nod[elnod[3]],x,y) &&
	 isnodonlhsof3pcurve (nod[elnod[3]],nod[elnod[7]],nod[elnod[0]],x,y)     )
      return true;
    break;
  }
  //case 413:{ break; }
  default:{
    print_err("wrong dimdegnne",__FILE__,__LINE__,__func__);
    exit (0);
    break;
  }}
  
  return false;
}


/**
   Function finds out element including "point".
   First of all, the elements in 'susel' are investigated.
   
   @param gt - gtopology, where the element is finding out
   @param nsusel - number suspected elements
   @param susel - array of suspected elements, size is gt->ne
   @param newsusel - empty array, size is gt->ne
   @param passedel - zero !!!!! array, size is gt->ne
   @param x,y - coordinates of point

   created  3.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
long whereispoint (gtopology *gt, long nsusel, long *susel, long *newsusel, long *passedel, double x, double y, long rep)
{
  long i,j,el,nnewsusel;
  
  /// point is in susel array
  if (rep==0)
    for (i=0;i<nsusel;i++)
      if ( !passedel[susel[i]]++ && (ispointinel (gt->gnodes,gt->gelements[susel[i]].nodes,x,y,gt->gelements[susel[i]].auxinf)) )
	return (susel[i]);
  
  nnewsusel = 0;
  for (i=0; i<nsusel; i++) {
    for (j=0; j<gt->nadjelel[susel[i]]; j++) {
      el = gt->adjelel[susel[i]][j];
      if (!passedel[el]){
	if (ispointinel (gt->gnodes,gt->gelements[el].nodes,x,y,gt->gelements[el].auxinf))
	  return (el);
	
	passedel[el] = 1;
	newsusel[nnewsusel++] = el;
      }
    }
  }
  
  if (!nnewsusel || rep==5) {
    print_err("Point does not lay in domain",__FILE__,__LINE__,__func__);
    return (-1);
  }
  
  return ( whereispoint (gt,nnewsusel,newsusel,susel,passedel,x,y,rep+1) );
}

double distance2 (double *a,double *b,long dim)
{
  if      (dim == 1)  return ((a[0]-b[0])*(a[0]-b[0]));
  else if (dim == 2)  return ((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]));
  else if (dim == 3)  return ((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]));
  else { print_err("Wrong dimension",__FILE__,__LINE__,__func__);  return (-1); }
  return 0;
}

/**
   Function finds out element including "point"
   or element closest to "point" in case the point doesnot lay in domain.
   
   @param dim - dimension of solved problem
   @param gt - gtopology, where the element is finding out
   @param x,y - coordinates of point
   
   created  11.12.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
long whereispoint_outpoint (gtopology *gt, long dim, double x, double y, double z)
{
  long i,cl,nadjel,*adjel;
  double aux, centr[3], point[3],dist=1e100;
  
  point[0]=x;
  point[1]=y;
  point[1]=z;
  
  //cl = gt->closest_node (dim,x,y);
  cl = closest_node (gt, dim, x,y,z);
  nadjel = gt->nadjelnod[cl];
  adjel = gt->adjelnod[cl];
  
  for (i=0;i<nadjel;i++)
    if ( ispointinel (gt->gnodes,gt->gelements[adjel[i]].nodes,x,y,gt->gelements[adjel[i]].auxinf) )
      return (adjel[i]);
  
  for (i=0;i<nadjel;i++){
    gt->gelements[adjel[i]].centroid (dim,gt->gnodes,centr);
    aux = distance2 (point,centr,dim);
    if (aux<dist){
      dist=aux;
      cl=i;
    }
  }
  
  return adjel[cl];
}

/**
   For every node (from new mesh) function finds out "parent" element (from old mesh),
   in which the "new" node lays.
   
   @param gt_old - gtopology of old mesh
   @param gt_new - gtopology of new mesh
   @param parentel - answer = allocated array of parent elements, size is mt_new->nn
   
   created  4.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void findout_parentel_nod (gtopology *gt_old, gtopology *gt_new, long dim, long *parentel)
{
  long i,j,k, n, nod, endheap;
  ivector oldsusel(gt_old->ne), newsusel(gt_old->ne), passedel(gt_old->ne), nodheap(gt_new->nn), nsusel(gt_new->nn);
  imatrix susel;
  
  if (dim==2)  allocm (gt_new->nn,18,susel);
  if (dim==3)  allocm (gt_new->nn,120,susel);
  
  fillv (-1, parentel, gt_new->nn);
  
  nodheap[0] = 0;
  endheap = 1;
  
  //nod = gt_old->closest_node (dim, gt_new->gnodes[nodheap[0]].x, gt_new->gnodes[nodheap[0]].y);
  nod = closest_node (gt_old, dim, gt_new->gnodes[nodheap[0]].x, gt_new->gnodes[nodheap[0]].y, gt_new->gnodes[nodheap[0]].z);
  
  nsusel[nodheap[0]] = gt_old->nadjelnod[nod];
  for (i=0; i<nsusel[nodheap[0]]; i++)
    susel[nodheap[0]][i] = gt_old->adjelnod[nod][i];
  
  for (i=0; i<gt_new->nn; i++) {
    nod = nodheap[i];
    
    nullv (passedel);
    copyv (susel[nod],oldsusel.a,nsusel[nod]);
    
    parentel[nod] = whereispoint (gt_old, nsusel[nod], oldsusel.a, newsusel.a, passedel.a, gt_new->gnodes[nod].x, gt_new->gnodes[nod].y, 0);
    if (parentel[nod]==-1)
      parentel[nod] = whereispoint_outpoint (gt_old, dim, gt_new->gnodes[nod].x, gt_new->gnodes[nod].y, gt_new->gnodes[nod].z);
    
    for (j=0;j<gt_new->nadjelnod[nod];j++)
      for (k=0;k<gt_new->gelements[gt_new->adjelnod[nod][j]].nne;k++){
	n = gt_new->gelements[gt_new->adjelnod[nod][j]].nodes[k];
	if (parentel[n] == -1){
	  parentel[n] = -2;
	  nodheap[endheap++] = n;
	}
	if (parentel[n] == -2)
	  susel[n][nsusel[n]++] = parentel[nod];
      }
    
  }
}

/**
   Function approximates 'values' from nodes into 'points'.
   All the point lays in 'element', which is determined by 'gt' and 'eid'.
   Values in points are returned by array pointval.
   
   @param npoints - number of points
   @param xx - array of x-coordinates of points, size is npoints
   @param yy - array of y-coordinates of points, size is npoints
   @param nval - number of values
   @param nodvalues - 1D array of values in nodes, size is nval*gt->nn, value position in array is
                      (val[0] on nod[0] ... val[0] on nod[gt->nn] ... val[nval] on nod[0] ... val[nval] on nod[gt->nn])
   @param pointval - answer = 2D array, size is npoints*nval
   @param pvmtrx - true == pointval is matrix; false == pointval is vector
   
   created  7.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void give_valuesinpoints (gtopology *gt, long eid, long npoints, double *xx, double *yy, long nval, const double *nodvalues, double **pointval, bool pvmtrx)
{
  long i,j,k,nne;
  double xi,eta,value;
  vector x, y, bf, nodval;
  ivector nodes;
  
  nne = gt->gelements[eid].nne;
  
  allocv (nne,x);
  allocv (nne,y);
  allocv (nne,bf);
  allocv (nne,nodes);
  allocv (nne,nodval);
  
  gt->give_nodes (eid,nodes);
  gt->give_node_coord2d (x,y,eid);
  
  for (i=0; i<npoints; i++) {
    
    switch (gt->gelements[eid].auxinf) {
    case 312:  nc_lin_3_2d  (     xx[i],yy[i],x.a,y.a, xi,eta);  bf_lin_3_2d  (bf.a,xi,eta);  break;
    case 622:  nc_quad_3_2d (1e-8,xx[i],yy[i],x.a,y.a, xi,eta);  bf_quad_3_2d (bf.a,xi,eta);  break;
    case 412:  nc_lin_4_2d  (1e-5,xx[i],yy[i],x.a,y.a, xi,eta);  bf_lin_4_2d  (bf.a,xi,eta);  break;
    case 822:  nc_quad_4_2d (1e-5,xx[i],yy[i],x.a,y.a, xi,eta);  bf_quad_4_2d (bf.a,xi,eta);  break;
    default: {
      print_err("unknown nnedegdim", __FILE__, __LINE__, __func__);  exit (1);
    }}
    
    for (j=0; j<nval; j++) {
      for (k=0; k<nne; k++)
	nodval[k] = nodvalues[nodes[k]+j*gt->nn];
      
      scprd (bf,nodval,value);
      
      if (pvmtrx)  pointval[i][j] = value;
      else         pointval[0][i+j*npoints] = value;
    }
  }
}

/**
   Function transforms values from nodes of old mesh into nodes of new mesh.
   
   @param po - problem of old mesh
   @param pn - problem of new mesh
   @param lcid - id of load case
   @param dim - dimension of solved problem
   @param ndofn - number of DOFs in node
   @param parentel - array of parent (old)elements of new mesh nodes, size is p_new->mt->nn
   @param spr
   
   created  7.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void transfvalues_nod (gtopology *gt_old, gtopology *gt_new, const double *r_old, double *r_new, long *parentel, long dim, long ndofn, bool spr)
{
  long i,j,k,nid,maxnchild;
  long *nchilds,**childnod;
  double **spcoord, *r_loc=NULL;
  vector x,y;
  
  if (dim != 2) { print_err("dim != 2", __FILE__, __LINE__, __func__);  exit (1); }
  
  /// CHILD NODES AT OLD ELEMENTS
  nchilds = new long [gt_old->ne];
  memset (nchilds, 0, gt_old->ne*sizeof(long));
  
  for (i=0; i<gt_new->nn; i++)
    nchilds[parentel[i]]++;
  
  childnod = new long* [gt_old->ne];
  for (i=0; i<gt_old->ne; i++)
    childnod[i] = new long [nchilds[i]];
  
  memset (nchilds,0,gt_old->ne*sizeof(long));
  for (i=0; i<gt_new->nn; i++)
    childnod[parentel[i]][nchilds[parentel[i]]++] = i;
  
  
  /// ALLOCATION
  if (spr) {
    spcoord = new double* [gt_old->ne];
    for (i=0; i<gt_old->ne; i++)  spcoord[i] = new double [nchilds[i]*dim];
  }
  else {
    maxnchild = 0;
    for (i=0; i<gt_old->ne; i++)  if (maxnchild<nchilds[i]) maxnchild = nchilds[i];
    allocv (maxnchild, x);
    allocv (maxnchild, y);
    r_loc = new double [maxnchild*ndofn];
  }
  
  // ********************
  /*
    double X,Y,f;
    
    for (i=0;i<po->mt->nn;i++){
    X = gt_old->gnodes[i].x;
    Y = gt_old->gnodes[i].y;
    f = -X*X*X*X/1000.0 - X*X*X/10.0 + X*X + 10*X -Y*Y*Y*Y/3.0 - Y*Y*Y/2.0 + 20*Y*Y + 30*Y - 212;
    
    r_old[i           ] = f;
    r_old[i+po->mt->nn] = f;
    }
    
    long hod,min;
    double sec = clock();
  */
  // ********************
  
  /// transformation from r_old to r_new
  if (spr) {
    { print_err("nezkontrolovano", __FILE__, __LINE__, __func__);  exit (1); }
    //for (i=0; i<gt_old->ne; i++) {
    //  for (j=0; j<nchilds[i]; j++) {
    //	spcoord[i][j*dim  ] = gt_new->gnodes[childnod[i][j]].x;
    //	spcoord[i][j*dim+1] = gt_new->gnodes[childnod[i][j]].y;
    //	if (dim==3) spcoord[i][j*dim+2] = gt_new->gnodes[childnod[i][j]].z;
    //  }
    //}
    //
    //least_square spr(dim,ndofn,gt_old,((Mp->adaptflag & 16) ? 16 : 0));
    //spr.L2_nod2sp (Out,nchilds,spcoord,r_old,r_new,'n');
  }
  else {
    for (i=0; i<gt_old->ne; i++) {
      for (j=0; j<nchilds[i]; j++) {
	x[j] = gt_new->gnodes[childnod[i][j]].x;
	y[j] = gt_new->gnodes[childnod[i][j]].y;
      }
      
      give_valuesinpoints (gt_old, i, nchilds[i], x.a, y.a, ndofn, r_old, &r_loc, false);
      
      // saving deformations from r_old to p_new->lsrs->lhs
      for (j=0; j<nchilds[i]; j++) {
	nid = childnod[i][j];
	for (k=0; k<ndofn; k++)
	  r_new[nid+k*gt_new->nn] = r_loc[j+k*nchilds[i]];
      }
    }
  }
  
  // ********************
  /*
    sec = (clock() - sec) / (double)CLOCKS_PER_SEC;
    hod = (long)sec/3600;  sec -= hod*3600;
    min = (long)sec/60;    sec -= min*60;
    fprintf (stdout,"\n -----------------------------------");
    fprintf (stdout,"\n Consumed time by TRANSF %2ld:%02ld:%05.2f",hod,min,sec);
    fprintf (stdout,"\n -----------------------------------\n");
    
    double a,ai,nordr,norr=0.0;
    vector dr(2*gt_new->nn),r(gt_new->nn);
    
    for (i=0;i<gt_old->ne;i++)
    for (j=0;j<nchilds[i];j++){
    nid = childnod[i][j];
    X = gt_new->gnodes[nid].x;
    Y = gt_new->gnodes[nid].y;
    r[nid] = -X*X*X*X/1000.0 - X*X*X/10.0 + X*X + 10*X -Y*Y*Y*Y/3.0 - Y*Y*Y/2.0 + 20*Y*Y + 30*Y - 212;
    
    dr[nid+gt_new->nn] = r[nid] - r_new[i][j*2+1];
    dr[nid           ] = r[nid] - r_new[i][j*2  ];
    }
    
    sizev (dr,nordr);
    fprintf (stdout,"\n\n nordr     =  %25.15f \n\n",nordr);
    fprintf (stdout,"\n\n nordr/nn1 =  %25.15f \n\n",nordr/sqrt(769));
    fprintf (stdout,"\n\n nordr/nn2 =  %25.15f \n\n",nordr/sqrt(24048));
    
    nordr = a = 0.0;
    for (i=0;i<gt_new->ne;i++){
    nordr += Pelt->error(i,dr,ai);
    a += ai;
    }
         
    nordr = sqrt(nordr/a);
    fprintf (stdout,"\n\n nordr(A)  =  %25.15f \n\n",nordr);
    
    exit (1);
    
    
    for (i=0;i<gt_new->nn;i++)  norr += r[i];
    sizev (r,norr); norr *= 2;
  */
  // ********************
  
  delete [] r_loc;
  delete [] nchilds;
  destrm (childnod,gt_old->ne);
  if (spr)  destrm (spcoord,gt_old->ne);
}



///
void adaptivityt :: statedata_backup (void)
{
  if (Tb->nlc != ntm) { print_err("Tb->nlc != dim", __FILE__, __LINE__, __func__);  exit (1); }
  
  // array of unknown values and derivatives at nodes
  delete [] r;      this->r     = new double [ntm * Gtt->nn];
  delete [] rdr;    this->rdr   = new double [ntm * Gtt->nn];
  
  long i, j, ii;
  for (i=0; i<Gtt->nn; i++)
    for (j=0; j<ntm; j++) {
      ii = Gtt->give_dof(i,j);
      if (ii>0) {
	this->r  [i+j*Gtt->nn] = Lsrst->  lhs[ii-1];
	this->rdr[i+j*Gtt->nn] = Lsrst->tdlhs[ii-1];
      }
      else {
	this->r  [i+j*Gtt->nn] = 0.0;
	this->rdr[i+j*Gtt->nn] = 0.0;
      }
    }
  
  ///
  delete    tctrl;  this->tctrl = new timecontr;                this->tctrl->be_copy_of (Tp->timecont);
  ;                                                             this->istep  =           Tp->istep;
}

///
void adaptivityt :: statedata_transfer (adaptivityt *Adat_old, gtopology *Gtt_old)
{
  /// transfer of r and rdr
  long *parentel = new long [Gtt->nn];
  findout_parentel_nod (Gtt_old, Gtt, dim, parentel);
  
  bool nod_spr = false;
  
  r   = new double [ntm * Gtt->nn];
  rdr = new double [ntm * Gtt->nn];
  transfvalues_nod (Gtt_old, Gtt, Adat_old->r,   r,   parentel, dim, ntm, nod_spr);
  transfvalues_nod (Gtt_old, Gtt, Adat_old->rdr, rdr, parentel, dim, ntm, nod_spr);
  
  delete [] parentel;
  
  
  /// transfer of tctrl
  if (tctrl)  { print_err("tctrl is not NULL",__FILE__,__LINE__,__func__); exit (0); }
  this->tctrl = Adat_old->tctrl;  Adat_old->tctrl = NULL;
  
  /// transfer of step
  istep = Adat_old->istep;
}

///
void adaptivityt :: statedata_restore (void)
{
  long i, j, ii;
  for (i=0; i<Gtt->nn; i++)
    for (j=0; j<ntm; j++) {
      ii = Gtt->give_dof(i,j);
      if (ii>0) {
	Lsrst->  lhs[ii-1] = this->r  [i+j*Gtt->nn];
	Lsrst->tdlhs[ii-1] = this->rdr[i+j*Gtt->nn];
      }
    }
  
  ///
  Tp->istep = this->istep;
}
