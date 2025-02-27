#include "meshtransfert.h"
#include "globalt.h"
#include "globmatt.h"
#include "gtopology.h"
#include "gelement.h"
#include "gnode.h"
#include "gmatrix.h"
#include "mathem.h"
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "least_square.h"
/**
   Function approximates values stored on integration points to nodes.
   Used interpolation method is 'spr_smoothing'.
   Values in nodes are returned by array nodvalues.
   Value position in array is (val on nod[0] ... val on nod[gt->nn])
   gt->nadjelem have to be allocated == to use function gt->adjacelem.
   
   @param ipvalues - values in integration points - array, dimension is total number of ip (tinp)

   @param nodvalues - empty(returned) array, dimension is number of nodes (Tt->nn)
   
   created  9.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
   modified by Tomas Krejci, 13.11.2003
 */
/*
void give_nodvals_ip (double *ipvalues,double *nodvalues)
{
  long i,j,k,nip,ipp,ntm,*nsp,dim;
  double **spcoord,**spvalue,**ipcoord;
  
  nsp = new long [Tt->ne];
  spcoord = new double* [Tt->ne];
  spvalue = new double* [Tt->ne];

  ntm = Tp->ntm;

  ipcoord = new double* [Tm->tnip];
  for (i=0;i<Tm->tnip;i++)
    ipcoord[i] = new double [3];
  
  allipcoord (ipcoord);
  
  for (i=0;i<Tt->ne;i++){
    nip = 0;
    for (j=0;j<ntm;j++)
      for (k=0;k<ntm;k++)
	nip += Tt->give_nip(i,j,k);
    
    nsp[i] = nip;
    dim = Tt->give_dimension (i);
    
    ipp = Tt->elements[i].ipp[0][0];
    spcoord[i] = new double [nip*dim];
    spvalue[i] = new double [nip];
    for (j=0;j<nip;j++){
      for (k=0;k<dim;k++)
	spcoord[i][j*dim+k] = ipcoord[ipp+j][k];
      spvalue[i][j] = ipvalues[ipp+j];
    }
  }
  
  least_square spr_est(dim,1,Gtt,1);
  
  spr_est.run (Outt,nsp,spcoord,spvalue,nodvalues);
  
  delete [] nsp;
  for (i=0;i<Tt->ne;i++){ delete [] spcoord[i]; delete [] spvalue[i]; }
  delete [] spcoord;
  delete [] spvalue;
  for (i=0;i<Tm->tnip;i++) delete [] ipcoord[i];
  delete [] ipcoord;
}
*/

/**
   function computes coordinates of all integration points in "global" mesh
   
   @param ipcoord - empty array, dimension is Tt->tnip x 3
   
   created  5.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
   modified by Tomas Krejci 13.11.2003
 */
/*
void allipcoord (double **ipcoord)
{
  long i,ri,ci,ntm,ipp;
  
  for (i=0;i<Tt->ne;i++){
    ntm = Tp->ntm;
    for (ri=0;ri<ntm;ri++){
      for (ci=0;ci<ntm;ci++){
	ipp = Tt->elements[i].ipp[ri][ci];
	
	switch (Tt->give_elem_type(i)){
	  //case barlint:{ Lbt->ipcoordblock (i,ri,ci,ipcoord+ipp); break; }
	  //case barlintax:{ Lbat->ipcoordblock (i,ri,ci,ipcoord+ipp); break; }
	  //case barquadt:{ Qbt->ipcoordblock (i,ri,ci,ipcoord+ipp); break; }
	  //case barquadtax:{ Qbat->ipcoordblock (i,ri,ci,ipcoord+ipp); break; }
	case trlint:{ Ltt->ipcoordblock (i,ri,ci,ipcoord+ipp); break; }
	case trlaxisym:{ Ltat->ipcoordblock (i,ri,ci,ipcoord+ipp); break; }
	case quadlint:{ Lqt->ipcoordblock (i,ri,ci,ipcoord+ipp); break; }
	case quadquadt:{ Qqt->ipcoordblock (i,ri,ci,ipcoord+ipp); break; }
	case quadquadtax:{ Qqat->ipcoordblock (i,ri,ci,ipcoord+ipp); break; }
	case quadlaxisym:{ Lqat->ipcoordblock (i,ri,ci,ipcoord+ipp); break; }
	default:{
	  fprintf (stderr,"\n\n unknown element type is required in function  allipcoord (%s, line %d).\n",__FILE__,__LINE__);
	}
	}
      }
    }
  }
}
*/


/**
   The function finds out node (from gt), which is the closest by required point with the given coordinates

   @param gt - pointer to the general topology class where the given mesh is stored
   @param x,y,z - coordinates of point
   
   created  30.1.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
 /*
long adjacnodet (gtopology *gt,double x,double y,double z)
{
  long i,a = 0;
  double d,min = 1e50;
  
  for(i=0;i<gt->ne;i++){
    d = (gt->gnodes[i].x - x)*(gt->gnodes[i].x - x) + (gt->gnodes[i].y - y)*(gt->gnodes[i].y - y) + (gt->gnodes[i].z - z)*(gt->gnodes[i].z - z);
    if (d<min){
      min=d;
      a=i;
    }
  }
  return (a);
}

 */
