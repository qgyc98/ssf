#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define EXTERN
#include "stackdef.h"
#include "globalg.h"

#include "cmlfile.h"
#include "optim_driver.h"

enum problem_type { optimization=1 } ;

#include "cr.h"
#include "csv.h"

//int main (int argc,char *argv[])
int main (int argc, const char *argv[])
{
  /*
  if ( argc<2 ) {
    printf( "usage: gefel <filename> \n" ) ;
    exit( 0 ) ;
  }

  time_t bt,et;
  long hod,min;
  double sec = clock();
  
  bt = time (NULL);
  
  //  program initiation and data reading
  long temp ;
  problem_type pt ;
  cmlfile *f ;
  
  f = new cmlfile( argv[1] ) ; 
  f->set_labels( 6 ) ;
  f->set_label_string( 0,"Problem" ) ;
  f->set_label_string( 1,"Algorithm" ) ;
  f->set_label_string( 2,"CallsLimit" ) ;
  f->set_label_string( 3,"Dim" ) ;
  f->set_label_string( 4,"Precision" ) ;
  f->set_label_string( 5,"ReturnToDomain" ) ;

  f->set_sections( 1 ) ;
  f->set_section_string( 0,"Domain" ) ;

  f->get_value(  0L, temp ) ;
  pt = ( problem_type ) temp ;
  fprintf (stdout,"\n Type of problem:         %d\n", pt );
  
  switch (pt)
    {
    case optimization:
      { 
	long calls_limit ;
	long dim ;
	algorithm_type at ;
	optim_driver *od ;

	f->get_value( 1, temp ) ;
	at = ( algorithm_type ) temp ;
	f->get_value( 2, calls_limit ) ;
	f->get_value( 3, dim ) ;
	fprintf (stdout," Type of algorithm   :    %d \n", at ) ;
	fprintf (stdout," Function calls limit:    %ld \n", calls_limit ) ;
	fprintf (stdout," Number of variables :    %ld \n", dim ) ;
	
	obj_funct F( dim ) ;

	long lines ;
	f->find_section( 0, lines ) ;
	if ( lines != dim ) {
	  fprintf (stderr,"\n\n Mismatch in definition of Domain") ;
	  fprintf (stderr,"\n main (file %s, line %d).\n",
		   __FILE__,__LINE__) ;
	  exit ( 1 ) ;
	}
	
	for ( long i=0 ; i<dim ; i++ ) {
	  f->get_section_value ( 0, i, 0, F.Domain[i][0] ) ;
	  f->get_section_value ( 0, i, 1, F.Domain[i][1] ) ;
	  f->get_section_value ( 0, i, 2, F.Domain[i][2] ) ;
	}
	f->end_reading_section() ;

	f->get_value( 4, F.precision ) ;
	f->get_value( 5, F.Return_to_domain ) ;

	fprintf (stdout," Precision of optimization is:    %f \n",\
		 F.precision ) ;
	fprintf (stdout," Return_to_domain is set to:    %ld \n",\
		 F.Return_to_domain ) ;

	od = new optim_driver( at, &F ) ;
	od->run( calls_limit ) ;
	break ; 
      }
    default:
      fprintf (stderr,"\n\n unknown type of problem in");
      fprintf (stderr,"\n main (file %s, line %d).\n",
	       __FILE__,__LINE__);
      exit ( 1 ) ;
    }

  f->close() ;
  delete f ;
  

  //  solution of the problem
  //  the most important function in the code

  
  et = time (NULL);
  fprintf (stdout,"\n\n\n Data about computation time \n");
  fprintf (stdout,"\n\n Total time of computation             %ld",et-bt);

  sec = (clock() - sec) / (double)CLOCKS_PER_SEC;
  hod = (long)sec/3600;  sec -= hod*3600;
  min = (long)sec/60;    sec -= min*60;
  fprintf (stdout,"\n ----------------------------------");
  fprintf (stdout,"\n Consumed time by GEFEL %2ld:%02ld:%05.2f",hod,min,sec);
  fprintf (stdout,"\n ----------------------------------\n");
  
  fprintf (stdout,"\n\n");  fprintf (stderr,"\n\n");  fprintf (stdout,"\n");
  */
  
  
  
  /*
    toto jsou testovaci casti od JK
  comprow *crm;
  compvect *cvi,*cvo;
  
  crm = new comprow;
  cvi = new compvect;
  cvo = new compvect;
  
  crm->n=5;
  crm->negm=15;
  crm->a=new double [15];
  crm->ci=new long [15];
  crm->adr=new long[6];
  
  crm->a[0]=1.0;
  crm->a[1]=2.0;
  crm->a[2]=3.0;
  crm->a[3]=4.0;
  crm->a[4]=5.0;
  crm->a[5]=6.0;
  crm->a[6]=7.0;
  crm->a[7]=1.0;
  crm->a[8]=2.0;
  crm->a[9]=3.0;
  crm->a[10]=4.0;
  crm->a[11]=5.0;
  crm->a[12]=6.0;
  crm->a[13]=7.0;
  crm->a[14]=8.0;
  
  crm->ci[0]=0;
  crm->ci[1]=1;
  crm->ci[2]=3;
  crm->ci[3]=4;
  crm->ci[4]=0;
  crm->ci[5]=2;
  crm->ci[6]=4;
  crm->ci[7]=0;
  crm->ci[8]=1;
  crm->ci[9]=2;
  crm->ci[10]=3;
  crm->ci[11]=4;
  crm->ci[12]=1;
  crm->ci[13]=2;
  crm->ci[14]=3;
  
  crm->adr[0]=0;
  crm->adr[1]=2;
  crm->adr[2]=4;
  crm->adr[3]=7;
  crm->adr[4]=12;
  crm->adr[5]=15;
  
  cvi->n=5;
  cvi->nz=2;

  cvi->a = new double [2];
  cvi->a[0]=3.0;
  cvi->a[1]=1.0;
  
  cvi->ind = new long [2];
  cvi->ind[0]=0;
  cvi->ind[1]=3;
  
  crm->crxcv_cv (cvi,cvo);
  
  long *li;
  long nsub=3;
  comprow *sm;
  
  sm = new comprow;
  li=new long[3];
  li[0]=1;
  li[1]=3;
  li[2]=4;

  crm->select_submatrix (li,nsub,sm);
  */
  
  long i,j,nn,ne,nne,ndofe;
  gelemtype get;
  gtopology *gtm;
  XFILE *in;
  FILE *out;

  if (argc < 2)
    return 0;
  
  out = fopen ("vysl.txt","w");
  in = xfopen(argv[1],"r");
  in->warning = 1;
  in->kwdmode = ignore_kwd;
  in->ignorecase = 1;


  gtm = new gtopology;

  xfscanf (in,"%k%ld","number_of_nodes",&nn);
  fprintf (stdout,"\n number of nodes  %ld",nn);
  gtm->alloc_nodes (nn);

  for (i=0;i<nn;i++){
    xfscanf (in,"%ld",&j);
    gtm->gnodes[j-1].read (in);
  }

  xfscanf (in,"%k%ld","number_of_elements",&ne);
  fprintf (stdout,"\n number of elements  %ld",ne);
  gtm->alloc_elements (ne);

  for (i=0;i<ne;i++){
    xfscanf (in,"%ld",&j);
    ndofe=24;
    nne=8;
    get = linhexa;
    gtm->gelements[i].read (in,nne,ndofe,get,nn);
  }

  gtm->adjacnodes (out);
  gtm->adjacelem (out);
  
  fclose (out);
  fprintf (stdout,"\n");
  return 0;
}
