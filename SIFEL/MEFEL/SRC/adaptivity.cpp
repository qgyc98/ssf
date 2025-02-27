#include "adaptivity.h"
#include "gadaptivity.h"
#include "element.h"
#include "intpoints.h"
#include "z2_smoothing.h"

#include "global.h"
#include "elemhead.h"
#include "elemswitch.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechbclc.h"
#include "lhsrhs.h"


/// CONSTRUCTOR
adaptivity :: adaptivity ()
{
  tad        = Mp->adaptivityflag;
  printflags = 0;
  corr       = 1.0;
  otherflags = -1;
  
  path = filename = suffix = NULL;
  ni = new char[1];  ni[0] = '\0';
  
  refined_ders_in_nodes = NULL;
  
  answer = 0;
}
/// DESTRUCTOR
adaptivity :: ~adaptivity ()
{
  delete [] ni;
  delete [] path;
  delete [] filename;
  delete [] suffix;
}

/// read values from input file
void adaptivity :: readinit (XFILE *in)
{
  xfscanf (in, "%k%ld %k%lf %k%lf", "printflag", &printflags, "accuracy", &adapt_accuracy, "correction", &corr);
  
  if (Mp->adaptivityflag<1 || Mp->adaptivityflag>2)
    { print_err("unknown type of adaptivity is required",__FILE__,__LINE__,__func__); exit (0); }
  
  filename_decomposition (in->fname, path, filename, suffix);
}

/// print values to input file
void adaptivity :: printinit (FILE *out)
{
  fprintf (out,"%ld %ld %.3lf %.3lf\n", tad, printflags, adapt_accuracy, corr);
}

/**
   This is main function of class adaptivity, it runs in three steps:
   1. preparing strain
   2. computing of refined derivative (strain or stress, determined by otherflag) by spr_smoothing or z2_smoothing
   3. computing of error of FEM solution
   (case 3 = computing of accurate error)
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
long adaptivity :: run (long of, int width, long nincr)
{
  long hod,min;
  double sec=0.0;
  
  if (Mespr==1){
    fprintf (stdout,"\n\n *********************************");
    fprintf (stdout,  "\n ***   Adaptivity is started   ***");
    
    sec = clock();
  }
  
  ///
  dim =   Mt->give_dimension (0);
  ord =   Mt->give_degree (0);
  ncomp = Mt->give_ncomp (0,0);
  nn =    Mt->nn;
  ne =    Mt->ne;
  
  //  initialization of auxinf variable
  for (long i=0; i<ne; i++) {
    if (ord != Mt->give_degree    (i))  { print_err("mismatch in order", __FILE__, __LINE__, __func__);  exit (1); }
    if (dim != Mt->give_dimension (i))  { print_err("mismatch in dimension", __FILE__, __LINE__, __func__);  exit (1); }
    Gtm->gelements[i].auxinf = 100*Mt->give_nne(i) + 10*ord + dim;
  }
  
  otherflags = of;
  answer = 0;
  
  /// strainy jsou pocitany standardne nekde vejs, kduz se pocitaji nejaky forces
  //  compute_ipstrains (0);
  
  // KONTROLOVAT TO POMOCI TECHTO PROMENNYCH   
  // Mp->strainstate=1;
  // Mp->stressstate=1;
  
  //if (otherflags & 8)  refined_ders_in_nodes = new double [ncomp*nn*2];
  //else                 refined_ders_in_nodes = new double [ncomp*nn];
  refined_ders_in_nodes = new vector[nn];
  for (long i=0; i<nn; i++)  allocv (ncomp, refined_ders_in_nodes[i]);
  
  prepare_ni (width, nincr);
  this->check_consistency ();
  
  switch (tad){
  case 1:{
    spr ();
    compute_error ();
    break;
  }
  case 2:{
    z2_smoothing z2_est(ncomp);
    z2_est.run (refined_ders_in_nodes);
    compute_error ();
    break;
  }
  //case 3:{
  //  spr ();
  //  compute_error_fine ();
  //  break;
  //}
  default: print_err("unknown type of error smoothing is required in function", __FILE__, __LINE__, __func__);  exit (1);
  }
  
  ///
  //this->print_addat_vtk ();
  
  // *** print of test file ***
  if (printflags & 16)
    print_test ();
  
  
  delete [] refined_ders_in_nodes;
  
  if (Mespr==1){
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
void adaptivity :: prepare_ni (
                               int width,
                               long // nincr
                              )
{
  ni = new char[1+width+1];
  
  int n = atol(suffix+1);
  
  if (n)  sprintf (ni,".%0*d", width, n);
  else    ni[0] = '\0';
  
  //
  if (strncmp(filename, "this.sifel", 10) == 0)
    filename[4] = '\0';
}


/**
   Function checks consistency of the problem.
   
   created  3.4.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void adaptivity :: check_consistency (void) const
{
  int te = Mt->give_elem_type (0);
  for (long i=1; i<ne; i++)
    if (Mt->give_elem_type(i) != te) {
      if (te == planeelementlt && Mt->give_elem_type(i) == planeelementlq) continue;
      if (te == planeelementlq && Mt->give_elem_type(i) == planeelementlt) continue;
      print_err("Uniform type of elements is required, element %ld differs", __FILE__, __LINE__, __func__, i+1);
      exit (1);
    }
  
  // tady checkovat ze prvni blok ipp ma vice ipp nez druhy
  switch (te){
  case planeelementlt:
  case planeelementlq:
    if (Pelt && (Pelt->nb!=1 || Pelt->nip[0][0]!=1 ||                                                                   Pelt->ncomp[0]!=3)                     ) { print_err("Pelt has nonstandard intpoints", __FILE__, __LINE__, __func__);  exit (1); }
    if (Pelq && (Pelq->nb!=2 || Pelq->nip[0][0]!=4 || Pelq->nip[1][1]!=1 || Pelq->nip[1][0]!=0 || Pelq->nip[0][1]!=0 || Pelq->ncomp[0]!=2 || Pelq->ncomp[1]!=1)) { print_err("Pelq has nonstandard intpoints", __FILE__, __LINE__, __func__);  exit (1); }
    break;
  case planeelementqt:
  case planeelementqq:
    if (Peqt && (Peqt->nb!=2 || Peqt->nip[0][0]!=3 || Peqt->nip[1][1]!=3 || Peqt->nip[1][0]!=0 || Peqt->nip[0][1]!=0 || Peqt->ncomp[0]!=2 || Peqt->ncomp[1]!=1)) { print_err("Peqt has nonstandard intpoints", __FILE__, __LINE__, __func__);  exit (1); }
    if (Peqq && (Peqq->nb!=2 || Peqq->nip[0][0]!=4 || Peqq->nip[1][1]!=4 || Peqq->nip[1][0]!=0 || Peqq->nip[0][1]!=0 || Peqq->ncomp[0]!=2 || Peqq->ncomp[1]!=1)) { print_err("Peqq has nonstandard intpoints", __FILE__, __LINE__, __func__);  exit (1); }
    break;
  case cctel:
    if (Cct  && (Cct ->nb!=2 || Cct ->nip[0][0]!=1 || Cct ->nip[1][1]!=1 || Cct ->nip[1][0]!=0 || Cct ->nip[0][1]!=0 || Cct ->ncomp[0]!=3 || Cct ->ncomp[1]!=2)) { print_err("Cct  has nonstandard intpoints", __FILE__, __LINE__, __func__);  exit (1); }
    break;
  case lineartet:
    if (Ltet && (Ltet->nb!=1 || Ltet->nip[0][0]!=1 ||                                                                   Ltet->ncomp[0]!=6 )                    ) { print_err("Ltet has nonstandard intpoints", __FILE__, __LINE__, __func__);  exit (1); }
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
void adaptivity :: spr (void)
{
  // This loop prepares array 'spder' for using in function `least_square::solve`
  // It allocates and fills up array 'spder' by derivatives(strain or stress) in main Gauss's points.
  // There is no choice now, strain is prefered.
  matrix *spder = new matrix[ne];
  
  int nip;
  long i, j, ipp;
  vector vect;
  vect.n = ncomp;
  for (i=0; i<ne; i++) {
    // nip = Mt->give_nip(i, 0, 0); // prvni blok
    // tato fce zatim neni, nahore se nip kontroluje natvrdo
    nip = 0;
    switch (Mt->elements[i].te) {
    case planeelementlt:   nip = 1;  break;
    case planeelementlq:   { Pelq->elchar  (i,spder[i]); break;}
    case planeelementqt:   { Peqt->elchar  (i,spder[i]); break;}
    case planeelementqq:   { Peqq->elchar  (i,spder[i]); break;}
    case lineartet:        { Ltet->elchar  (i,spder[i]); break;}
    default: print_err("unknown element type", __FILE__, __LINE__, __func__);  exit (0);
    }
    
    if (nip) {
      allocm(nip, ncomp, spder[i]);
      ipp = Mt->elements[i].ipp[0][0];
      
      for (j=0; j<nip; j++) {
	vect.a = spder[i].a + j*ncomp;
	
	if (!(otherflags & 2))  Mm->givestrain (0, ipp, vect);
	else                    Mm->givestress (0,ipp,vect);
	ipp++;
      }
    }
  }
  vect.a = NULL;
  
  //spr.spr_default (Out,spder,refined_ders_in_nodes);
  patch_averaging spr (Gtm, dim, ncomp, 0);
  
  spr.solve (Out, spder, refined_ders_in_nodes);
  
  delete [] spder;
}

/**
   Function computes error of FEM solution from difference between
   apriximate derivative fied(from FEM solution) and much more exact derivative fied(from array refined_ders_in_nodes).
   Considering 'printflags' printing of any files and information is performed.
   
   created  3.6.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void adaptivity :: compute_error (void)
{
  long i;
  double e2     ,u2     ,volume   ,corr1,corr2;
  vector ei2(ne),ui2(ne),sizel(ne);
  
  e2 = u2 = volume = 0.0;
  
  for (i=0; i<ne; i++) {
    
    switch (Mt->give_elem_type (i)) {
    case planeelementlt:  volume += Pelt->compute_error (i, ei2[i], ui2[i], sizel[i], refined_ders_in_nodes,otherflags);  corr1 = 1.00; corr2 = 1.11;   break;
    case planeelementlq:  volume += Pelq->compute_error (i, ei2[i], ui2[i], sizel[i], refined_ders_in_nodes);  corr1 = 1.00; corr2 = 1.20;   break;
    case planeelementqt:  volume += Peqt->compute_error (i, ei2[i], ui2[i], sizel[i], refined_ders_in_nodes);  corr1 = 1.07; corr2 = 1.55;   break;
    case planeelementqq:  volume += Peqq->compute_error (i, ei2[i], ui2[i], sizel[i], refined_ders_in_nodes);  corr1 = 1.00; corr2 = 1.60;   break;
    case lineartet:       volume += Ltet->compute_error (i, ei2[i], ui2[i], sizel[i], refined_ders_in_nodes);  corr1 = 1.00; corr2 = 1.00;   break;
    default:  print_err("wrong mesh - element %ld", __FILE__, __LINE__, __func__, i+1);  exit (1);
    }
    
    if (tad == 1) ei2[i] *= corr1*corr1;
    if (tad == 2) ei2[i] *= corr2*corr2;
    
    e2 += ei2[i];
    u2 += ui2[i];
  }
  u2 += e2;
  
  //
  compute_refsizel_lin (sizel.a, ei2.a, e2, u2);
  
  /// percentual global error
  double er = sqrt(e2/u2);
  if (er > adapt_accuracy)  answer = 1;
  
  if (Mespr==1)
    fprintf (stdout, "\n global error = %.1lf%% (sqrt(e2/u2) == sqrt(%.3e / %.3e)) at volume = %f", 100*er, e2, u2, volume);
  
  /// *** print of background mesh file ***
  if (printflags & 2  &&  answer){
    char name[255];
    sprintf (name,"%s%s.T3d.bgm%s", path, filename, suffix);
    print_valel (Out, Gtm, path, name, "chyba", refsizel.a, 'e');
  }
  
  ///
  allocv (ne, elem_error_pct);
  for (i=0; i<ne; i++)
    elem_error_pct[i] = 100.0 * sqrt( ei2[i] / (ui2[i]+ei2[i]) );
  
}


/**
   Function computes new=required characteristic size of elements.
   
   @param sizel - array of old characteristic sizes of element
   @param
   
   created  3.6.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void adaptivity :: compute_refsizel_lin (const double *sizel, const double *ei2, double /*e2*/, double u2)
{
  long i;
  double em2,specacc;
  vector epsil;
  
  allocv (ne, epsil);
  
  specacc = this->adapt_accuracy;
  em2 = specacc*specacc * u2/ne;
  
  //corr = 0.8;  //<0.7;1.0>
  
  allocv (ne, refsizel);
  allocv (ne, refsizelrel);
  
  for (i=0; i<ne; i++) {
    epsil[i] = sqrt(ei2[i] / em2);
    if (epsil[i]<0.05) epsil[i] = 0.05;
    
    refsizelrel[i] = 1.0      / pow(epsil[i] , 1.0/(double)ord * corr);
    refsizel[i]    = sizel[i] / pow(epsil[i] , 1.0/(double)ord * corr);
  }
  
  // ***********************************************
  //char file[255];
  //sprintf (file,"%sepsil.dx%s",Mp->path,ni);
  //print_valel (Gtm,Mp,file,"epsil",epsil.a,'d');
  // ***********************************************
}



/**

   created  3.6.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void adaptivity :: print_test (void) const
{
  long i,j,ipp;
  char f_name[255];
  FILE *f_test;
  
  sprintf (f_name,"%s%s.sifel.test",path,filename);
  f_test = fopen (f_name,"a");            if (f_test==NULL) fprintf (stderr,"\n .test file has not been opened.\n\n");
  
  
  // ***************************************************************************  
  fprintf (f_test,"\n\n LEFT HAND SIDE \n\n");
  
  for (i=0;i<Mb->nlc*Ndofm;i++){
    fprintf (f_test," %20.10f ",Lsrs->lhs[i]);
    if ((i+1)%7==0)  fprintf (f_test,"\n");
  }
  
  
  // ***************************************************************************  
  fprintf (f_test,"\n\n STRAINS AND STRESSES ON INTEGRATION POINTS \n\n");
  
  for (i=0;i<ne;i++){
    ipp = Mt->elements[i].ipp[0][0];
    
    switch ((long)Mt->elements[i].te){
    case planeelementlq:{
      for (j=0;j<4;j++){
	fprintf (f_test," %20.10f "  ,Mm->ip[ipp+j].strain[0]);
	fprintf (f_test," %20.10f "  ,Mm->ip[ipp+j].strain[1]);
	fprintf (f_test," %20.10f \n",Mm->ip[ipp+4].strain[2]);
      }
      break;
    }
    case planeelementqq:{
      for (j=0;j<4;j++){
	fprintf (f_test," %20.10f "  ,Mm->ip[ipp+j  ].strain[0]);
	fprintf (f_test," %20.10f "  ,Mm->ip[ipp+j  ].strain[1]);
	fprintf (f_test," %20.10f \n",Mm->ip[ipp+j+4].strain[2]);
	}
      break;
    }
    case planeelementlt:{
      fprintf (f_test," %20.10f "  ,Mm->ip[ipp].strain[0]);
      fprintf (f_test," %20.10f "  ,Mm->ip[ipp].strain[1]);
      fprintf (f_test," %20.10f \n",Mm->ip[ipp].strain[2]);
      break;
    }
    case planeelementqt:
    case planeelementsubqt:{
      for (j=0;j<3;j++){
	fprintf (f_test," %20.10f "  ,Mm->ip[ipp+j  ].strain[0]);
	fprintf (f_test," %20.10f "  ,Mm->ip[ipp+j  ].strain[1]);
	fprintf (f_test," %20.10f \n",Mm->ip[ipp+j+3].strain[2]);
      }
      break;
    }
    case lineartet:{
      for (j=0;j<6;j++)
	fprintf (f_test," %20.10f ",Mm->ip[ipp].strain[j]);
      
      fprintf (f_test,"\n");
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown element type is required in function");
      fprintf (stderr," adaptivity::print_test  (%s, line %d)",__FILE__,__LINE__);
    }
    }
    
  }
  
  
  // ***************************************************************************  
  fprintf (f_test,"\n\n REFINED_DERS_IN_NODES \n\n");
  
  long l=0;
  for (i=0; i<ncomp; i++)
    for (j=0; j<nn; j++) {
      fprintf (f_test," %20.10f ", refined_ders_in_nodes[j].a[i]);
      if ((l+1)%7==0)  fprintf (f_test,"\n");
      l++;
    }
  
  //if (otherflags & 8){
  //  fprintf (f_test,"\n");
  //  
  //  for (i=0;i<ncomp*nn;i++){
  //    fprintf (f_test," %20.10f ", refined_ders_in_nodes[i + ncomp*nn]);
  //    if ((i+1)%7==0)  fprintf (f_test,"\n");
  //  }
  //}
  
  // ***************************************************************************  
  fprintf (f_test,"\n\n REFSIZEL \n\n");
  
  for (i=0;i<Mt->ne;i++){
    fprintf (f_test," %20.10f ",refsizel[i]);
    if ((i+1)%7==0)  fprintf (f_test,"\n");
  }
  
  // ***************************************************************************  
  if (Mm->ip[0].eqother!=NULL){
    fprintf (f_test,"\n\n OTHER \n\n");
    
    for (i=0;i<Mm->tnip;i++){
      for (j=0;j<6;j++)
	fprintf (f_test," %20.10f "  ,Mm->ip[i].eqother[j]);
      
      fprintf (f_test,"\n");
    }
  }
  
  fclose (f_test);
}
