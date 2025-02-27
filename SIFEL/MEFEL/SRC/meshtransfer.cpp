#include "meshtransfer.h"

#include "global.h"
#include "globmat.h"
#include "intpoints.h"
#include "loadcase.h"
#include "element.h"
#include "plelemlt.h"
#include "plelemqt.h"
#include "plelemlq.h"
#include "plelemqq.h"
#include "nssolver.h"
#include "stochdriver.h"
#include "mefelinit.h"

#include "gtopology.h"
#include "gelement.h"
#include "gnode.h"
#include "gmatrix.h"
#include "mathem.h"
#include "problem.h"
//#include "least_square.h"

#include <time.h>
#include <math.h>
#include <stdlib.h>

////////////////////       /* termitovo */       ////////////////////////////////////

//  /*
//  vse je delano pro 2D planestress prvky s ndofn=2
//  pro jine to neni kontrolovano !!!!!!!!!!!!!!!!!!
//  prvek plelemqq je predpokladan s rovnymi hranami(subparametricky), neb v fci 'nc_quad_4_2d' je zatim odkaz jen na 'nc_lin_4_2d'
//  */
//  /*
//  ze site na sit se prenaseji posuny a prvnich ncompstr hodnot z pole other , ne prvnich ncompother neb gamma a hardening nepotrebuju
//  */
//  /*
//  fce nezavisle na jine globalni promenne:
//  give_ndofn             ok
//  give_dof               ok
//  give_nip               ok
//  give_nne               ok
//  gt->give_nodes         ok
//  gt->give_node_coord2d  ok
//  */
//  /**
//   */
//  void newmeshread (const char *filename,long lcid)
//  {
//    long ip_dir,nod_spr;
//    ip_dir = 0;
//    nod_spr = 1;
//    
//    if ((Mp->adaptflag & 16))
//      nod_spr = 0;
//    
//    
//    // **********************************************************************
//    // begin of step 1  -  aproximation `values` into nodes
//    // **********************************************************************
//    
//    long ndofn,dim,nipcomp[3];
//    double *nodvals_ip;
//    
//    dim = 2;
//    ndofn = 2;
//    
//    nipcomp[0] = 0;   // number of transferd components of array ip[]->strain
//    nipcomp[1] = 0;   // number of transferd components of array ip[]->stress
//    nipcomp[2] = 5;   // number of transferd components of array ip[]->eqother        !!!  bez hardeningu
//    
//    if (!ip_dir)
//      give_nodvals_ip (dim,nipcomp,nodvals_ip);
//    else
//      nodvals_ip = NULL;
//    
//    // **********************************************************************
//    // begin of step 2  -  loading of new mesh
//    // **********************************************************************
//    
//    problem *p_old = new problem;
//    problem *p_new = new problem;
//    
//    p_old->globinic ();
//    
//    Mp = NULL; Gtm = NULL; Mt = NULL; Mm = NULL; Mc = NULL; Mb = NULL; Lsrs = NULL; Smat = NULL;
//    fclose (Out);
//    
//    stochdriver stochd;
//    const char *argv[2];
//    argv[0] = "";
//    argv[1] = filename;
//    mefel_init (2, argv, &stochd);
//    
//    p_new->globinic ();
//    
//    // **********************************************************************
//    // begin of step 3  -  konverse
//    // **********************************************************************
//    
//    long *parentel_ip,*parentel_nod;
//    double **ipcoord_new;
//    
//    parentel_ip = new long [p_new->mm->tnip];
//    parentel_nod = new long [p_new->mt->nn];
//    reallocm (p_new->mm->tnip,3,ipcoord_new);
//    
//    
//    findout_parentel_ip (p_old->gt,Gtm,Mm,Mt,ipcoord_new,parentel_ip);
//    findout_parentel_nod (p_old->gt,Gtm,parentel_nod,dim);
//    
//    //if (Mespr==1)   fprintf (stdout,"\n ***  Transformation of meshes  ***\n");
//    
//    if (!ip_dir)  transfvalues_ip_indirect (p_old->gt,Mm,                 ipcoord_new,parentel_ip,nipcomp,nodvals_ip );
//    else          transfvalues_ip_direct   (p_old->gt,Mm,(const double **)ipcoord_new,parentel_ip,nipcomp,p_old->mm,dim);
//    
//    if (!nod_spr) transfvalues_nod (p_old,p_new,lcid,dim,ndofn,parentel_nod,'n');
//    else          transfvalues_nod (p_old,p_new,lcid,dim,ndofn,parentel_nod,'y');
//    
//    
//    delete [] parentel_ip;
//    delete [] parentel_nod;
//    
//    // **********************************************************************
//    // end of step 3  -  konverse
//    // **********************************************************************
//    
//    p_old->dealoc ();
//    p_new->deinic ();
//    
//    // **********************************************************************
//    // end of step 2  -  loading of new mesh
//    // **********************************************************************
//    
//    if (nodvals_ip!=NULL)  delete [] nodvals_ip;
//    
//    // **********************************************************************
//    // end of step 1  -  aproximation `values` into nodes
//    // **********************************************************************
//    
//    
//    mefel_right_hand_side (lcid,Lsrs->rhs);
//  }
//  
//  /**
//     Function approximates values stored on integration points (in arrays `strain` `stress` `eqother`) to nodes - for "global" mesh.
//     Used interpolation method is 'spr_smoothing'.
//     Values in nodes are returned by array nodvalues.
//     Value position in array is (val[0] on nod[0] ... val[0] on nod[gt->nn] ... val[nval] on nod[0] ... val[nval] on nod[gt->nn])
//     Necessary precondition is Mm->ip[i].ncompstr and Mm->ip[i].ncompother are constant all over the domain !!!!
//     
//     @param dim - dimension of solved problem
//     @param nipcomp - number of members (of arrays Mm->ip[].strain , Mm->ip[].stress and Mm->ip[].eqother), which are approximated;
//                      nipcomp[0] <0;Mm->ip[].ncompstr>, nipcomp[1] <0;Mm->ip[].ncompstr>, nipcomp[2] <0;Mm->ip[].ncompother>
//  		    if nipcomp = -1 => all members
//     @param nodvalues - answer ; nonalocated
//     
//     created  9.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  long give_nodvals_ip (long dim,long nipcomp[],double *&nodvalues)
//  {
//    long i,j,k,nip,ipp,nvals,*nsp;
//    double **spcoord,**spvalue,**ipcoord;
//    
//    if (nipcomp[0]==-1 || nipcomp[0] > Mm->ip[0].ncompstr  ) nipcomp[0] = Mm->ip[0].ncompstr;
//    if (nipcomp[1]==-1 || nipcomp[1] > Mm->ip[0].ncompstr  ) nipcomp[1] = Mm->ip[0].ncompstr;
//    if (nipcomp[2]==-1 || nipcomp[2] > Mm->ip[0].ncompother) nipcomp[2] = Mm->ip[0].ncompother;
//    
//    nvals = nipcomp[0] + nipcomp[1] + nipcomp[2];
//    
//    if (nvals==0){
//      nodvalues = NULL;
//      return 0;
//    }
//    else
//      nodvalues = new double [nvals*Mt->nn];
//    
//    nsp = new long [Mt->ne];
//    spcoord = new double* [Mt->ne];
//    spvalue = new double* [Mt->ne];
//    ipcoord = new double* [Mm->tnip];
//    for (i=0;i<Mm->tnip;i++)
//      ipcoord[i] = new double [3];
//    
//    allipcoord (ipcoord);
//    
//    for (i=0;i<Mt->ne;i++){
//      nip = 0;
//      for (j=0;j<Mt->elements[i].nb;j++)
//        for (k=0;k<Mt->elements[i].nb;k++)
//  	nip += Mt->elements[i].nip[j][k];
//      
//      nsp[i] = nip;
//      
//      ipp = Mt->elements[i].ipp[0][0];
//      spcoord[i] = new double [nip*dim];
//      spvalue[i] = new double [nip*nvals];
//      for (j=0;j<nip;j++){
//        for (k=0;k<dim;k++)
//  	spcoord[i][j*dim+k] = ipcoord[ipp+j][k];
//        
//        for (k=0;k<nipcomp[0];k++)
//  	spvalue[i][j*nvals + k                          ] = Mm->ip[ipp+j].strain[k];
//        
//        for (k=0;k<nipcomp[1];k++)
//  	spvalue[i][j*nvals + k + nipcomp[0]             ] = Mm->ip[ipp+j].stress[k];
//        
//        for (k=0;k<nipcomp[2];k++)
//  	spvalue[i][j*nvals + k + nipcomp[0] + nipcomp[1]] = Mm->ip[ipp+j].eqother[k];
//        
//      }
//    }
//    
//    least_square spr(dim,nvals,Gtm,0+((Mp->adaptflag & 16) ? 16 : 0));   // co to je 1 ??
//    
//    spr.spr_optional (Out,nsp,spcoord,spvalue,nodvalues,1);
//    
//    delete [] nsp;
//    for (i=0;i<Mt->ne;i++){ delete [] spcoord[i]; delete [] spvalue[i]; }
//    delete [] spcoord;
//    delete [] spvalue;
//    for (i=0;i<Mm->tnip;i++) delete [] ipcoord[i];
//    delete [] ipcoord;
//    
//    return nvals;
//  }
//  
//  /**
//     Function returns array of deformations in ALL nodes of problem 'p'.
//     Deformation position in answer array is (r[node_1][x] ... r[node_nn][x] , r[node_1][y] ... r[node_nn][y])
//     
//     @param lcid  - id of load case
//     @param p     - pointr on problem
//     @param rfull - answer - allocated 1D array, size is p->mt->nn * p->mm->ip[0].ndofn
//     
//     created  9.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  void give_rfull (long lcid,problem *p,double *rfull)
//  {
//    long i,j,ndofn;
//    double *r;
//  
//    ndofn = p->mt->give_ndofn (0);
//    r = new double [ndofn];
//  
//    for (i=0;i<p->mt->nn;i++){
//      noddispl (lcid,i,p,r);
//      for (j=0;j<ndofn;j++)
//        rfull[i+j*p->mt->nn] = r[j];
//    }
//    
//    delete [] r;
//  }
//  
//  /**
//     function computes coordinates of all integration points in "global" mesh
//     
//     @param ipcoord - empty 1D array, size is Mt->tnip x 3
//     
//     created  5.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  void allipcoord (double **ipcoord)
//  {
//    long i,ri,ci,nb,ipp;
//    
//    for (i=0;i<Mt->ne;i++){
//      nb = Mt->elements[i].nb;
//      for (ri=0;ri<nb;ri++){
//        for (ci=0;ci<nb;ci++){
//  	ipp = Mt->elements[i].ipp[ri][ci];
//  	
//  	if (Mt->elements[i].nip[ri][ci]){
//  	  switch (Mt->give_elem_type(i)){
//  	  case planeelementlt:{ Pelt->ipcoordblock (i,ri,ci,ipcoord+ipp); break; }
//  	  case planeelementqt:{ Peqt->ipcoordblock (i,ri,ci,ipcoord+ipp); break; }
//  	  case planeelementlq:{ Pelq->ipcoordblock (i,ri,ci,ipcoord+ipp); break; }
//  	  case planeelementqq:{ Peqq->ipcoordblock (i,ri,ci,ipcoord+ipp); break; }
//  	  default:{
//  	    fprintf (stderr,"\n\n unknown element type is required in function  allipcoord (%s, line %d).\n",__FILE__,__LINE__);
//  	  }
//  	  }
//  	}
//        }
//      }
//    }
//  }
//  
//  /**
//     Function alocates required arrays for function 'locintpoints',
//     which finds out (for every integration point from new mesh) "parent" element (from old mesh),
//     in which the integration point lays.
//     Parent elements are returned in array "parentel_ip".
//     
//     @param gt_old - gtopology of old mesh
//     @param gt_new - gtopology of new mesh
//     @param mm_new - mechmat of new mesh
//     @param mt_new - mechtop of new mesh
//     @param ipcoord - array of coordinates of new mesh int. points, size is mm_new->tnip x dim
//     @param parentel_ip - answer = 1D array, size is mt_new->tnip
//     
//     created  3.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  void findout_parentel_ip (gtopology *gt_old,gtopology *gt_new,mechmat *mm_new,mechtop *mt_new,double **ipcoord,long *parentel_ip)
//  {
//    long *passedel1;  passedel1 = new long [gt_old->ne];
//    long *passedel2;  passedel2 = new long [gt_new->ne];   memset (passedel2,0,gt_new->ne*sizeof(long));
//    long *elheap;     elheap    = new long [gt_new->ne];   elheap[0] = 0;
//    long *newelheap;  newelheap = new long [gt_new->ne];
//    long *susel;      susel     = new long [gt_old->ne];
//    long *newsusel;   newsusel  = new long [gt_old->ne];
//    
//    memset (parentel_ip,0,mm_new->tnip*sizeof(long));
//    allipcoord (ipcoord);
//    
//    locintpoints (gt_old,gt_new,mt_new,ipcoord,passedel1,passedel2,1,elheap,newelheap,susel,newsusel,parentel_ip);
//    
//    delete [] passedel1;
//    delete [] passedel2;
//    delete [] elheap;
//    delete [] newelheap;
//    delete [] susel;
//    delete [] newsusel;
//  }
//  
//  /**
//     For every integration point (from new mesh) function finds out "parent" element (from old mesh),
//     in which the int. point lays.
//     Parent elements are returned in array "parentel".
//  
//     @param gt_old - gtopology of old mesh
//     @param gt_new - gtopology of new mesh
//     @param mt_new - mechtop of new mesh
//     @param ipcoord - array of coordinates of new mesh int. points, size is mm_new->tnip x dim
//     @param passedel1 - empty array, size is gt_old->ne
//     @param passedel2 - zero array, size is gt_new->ne
//     @param nelheap = 1
//     @param elheap - array of size gt_new->ne, elheap[0] = 0
//     @param newelheap - empty array, size is gt_new->ne
//     @param susel - empty array, size is gt_old->ne
//     @param newsusel - answer = 1D array, size is gt_old->ne
//     @param parentel - answer = 1D array, size is mm_new->tnip
//     
//     created  3.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  void locintpoints (gtopology *gt_old,gtopology *gt_new,mechtop *mt_new,double **ipcoord,long *passedel1,long *passedel2,long nelheap,long *elheap,long *newelheap,long *susel,long *newsusel,long *parentel)
//  {
//    long i,j,k,l,el,nip,ipp,nnewelheap,nsusel;
//    
//    nnewelheap = 0;
//    for (i=0;i<nelheap;i++){
//      nsusel = mt_new->give_nip(elheap[i],0,0);
//      ipp = mt_new->elements[elheap[i]].ipp[0][0];
//      for (j=0;j<nsusel;j++)
//        susel[j] = parentel[ipp++];
//      
//      for (j=0;j<gt_new->nadjelnod[elheap[i]];j++){
//        el = gt_new->adjelel[elheap[i]][j];
//        
//        if (!passedel2[el]){
//  	nip = 0;
//  	ipp = mt_new->elements[el].nb;
//  	for (k=0;k<ipp;k++)
//  	  for (l=0;l<ipp;l++)
//  	    nip += mt_new->elements[el].nip[k][l];
//  	
//  	ipp = mt_new->elements[el].ipp[0][0];
//  	for (k=0;k<nip;k++){
//  	  for (l=0;l<gt_old->ne;l++)
//  	    passedel1[l] = 0;
//  	  
//  	  parentel[ipp+k] = whereispoint (gt_old,nsusel,susel,newsusel,passedel1,ipcoord[ipp+k][0],ipcoord[ipp+k][1],'f');
//  	}
//  	
//  	passedel2[el] = 1;
//  	newelheap[nnewelheap++] = el;
//        }
//      }
//    }
//    
//    if (!nnewelheap)
//      return ;
//    
//    locintpoints (gt_old,gt_new,mt_new,ipcoord,passedel1,passedel2,nnewelheap,newelheap,elheap,susel,newsusel,parentel);
//    return ;
//  }
//  
//  /**
//     Function runs function 'locnodes'.
//     
//     @param gt_old - gtopology of old mesh
//     @param gt_new - gtopology of new mesh
//     @param parentel_nod - answer = array of parent elements, size is mt_new->nn
//     
//     created  6.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  void findout_parentel_nod (gtopology *gt_old,gtopology *gt_new,long *parentel_nod,long dim)
//  {
//    locnodes (gt_old,gt_new,parentel_nod,dim);
//  }
//  
//  /**
//     For every node (from new mesh) function finds out "parent" element (from old mesh),
//     in which the "new" node lays.
//     
//     @param gt_old - gtopology of old mesh
//     @param gt_new - gtopology of new mesh
//     @param parentel - answer = allocated array of parent elements, size is mt_new->nn
//     
//     created  4.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//  */
//  void locnodes (gtopology *gt_old,gtopology *gt_new,long *parentel,long dim)
//  {
//    long i,j,k,n,nod,endheap;
//    ivector oldsusel(gt_old->ne),newsusel(gt_old->ne),passedel(gt_old->ne),nodheap(gt_new->nn),nsusel(gt_new->nn);
//    imatrix susel;
//    
//    if (dim==2)  reallocm (gt_new->nn,18,susel);
//    if (dim==3)  reallocm (gt_new->nn,120,susel);
//    
//    fillv (-1,parentel,gt_new->nn);
//    
//    endheap = 1;
//    nsusel[0] = 1;
//    susel[0][0] = 0;
//    nodheap[0] = 0;
//    
//    for (i=0;i<gt_new->nn;i++){
//      nod = nodheap[i];
//      
//      nullv (passedel);
//      copyv (susel[nod],oldsusel.a,nsusel[nod]);
//      
//      parentel[nod] = whereispoint (gt_old,nsusel[nod],oldsusel.a,newsusel.a,passedel.a,gt_new->gnodes[nod].x,gt_new->gnodes[nod].y,'f');
//      
//      for (j=0;j<gt_new->nadjelnod[nod];j++)
//        for (k=0;k<gt_new->gelements[gt_new->adjelnod[nod][j]].nne;k++){
//  	n = gt_new->gelements[gt_new->adjelnod[nod][j]].nodes[k];
//  	if (parentel[n] == -1){
//  	  parentel[n] = -2;
//  	  nodheap[endheap++] = n;
//  	}
//  	if (parentel[n] == -2)
//  	  susel[n][nsusel[n]++] = parentel[nod];
//        }
//      
//    }
//  }
//  
//  /**
//     Function finds out element including "point".
//     First of all, the elements in "susel" are investigated.
//     
//     @param gt - gtopology, where the element is finding out
//     @param nsusel - number suspected elements
//     @param susel - array of suspected elements, size is gt->ne
//     @param newsusel - empty array, size is gt->ne
//     @param passedel - zero !!!!! array, size is gt->ne
//     @param x,y - coordinates of point
//  
//     created  3.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  long whereispoint (gtopology *gt,long nsusel,long *susel,long *newsusel,long *passedel,double x,double y,char flag)
//  {
//    long i,j,el,nnewsusel;
//    
//    if (flag=='f')
//      for (i=0;i<nsusel;i++)
//        if ( !passedel[susel[i]]++ && (ispointinel (gt->gnodes,gt->gelements[susel[i]].nodes,x,y,gt->gelements[susel[i]].auxinf)) )
//  	return (susel[i]);
//    
//    nnewsusel = 0;
//    for (i=0;i<nsusel;i++){
//      for (j=0;j<gt->nadjelnod[susel[i]];j++){
//        el = gt->adjelel[susel[i]][j];
//        if (!passedel[el]){
//  	if (ispointinel (gt->gnodes,gt->gelements[el].nodes,x,y,gt->gelements[el].auxinf))
//  	  return (el);
//  	
//  	passedel[el] = 1;
//  	newsusel[nnewsusel++] = el;
//        }
//      }
//    }
//    
//    if (!nnewsusel){
//      fprintf (stderr,"\n\n!!!!!!!!!  Point does not lay in domain (%s, line %d)\n\n",__FILE__,__LINE__);
//      return (-1);
//    }
//    
//    return ( whereispoint (gt,nnewsusel,newsusel,susel,passedel,x,y,'n') );
//  }
//  
//  /**
//     Function finds out, whether 'node' lays in 'element'.
//     Node is defined by its coordinates - x,y.
//     Element is defined by its node numbers - 'elnod' and dimdegnne - 'ndd'
//     
//     @param nod - array of all gnodes
//     @param elnod - array of node numbers of element
//     @param x,y - coordinates of node
//     @param ndd - nnedimdeg of element
//     
//     created  3.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  long ispointinel (gnode *nod,long *elnod,double x,double y,long ndd)
//  {
//    switch (ndd){
//    case 312:{
//      if ( isnodonlhsofline (nod[elnod[0]],nod[elnod[1]],x,y) &&
//  	 isnodonlhsofline (nod[elnod[1]],nod[elnod[2]],x,y) &&
//  	 isnodonlhsofline (nod[elnod[2]],nod[elnod[0]],x,y)     )
//        return (1);
//      break;
//    }
//    case 412:{
//      if ( isnodonlhsofline (nod[elnod[0]],nod[elnod[1]],x,y) &&
//  	 isnodonlhsofline (nod[elnod[1]],nod[elnod[2]],x,y) &&
//  	 isnodonlhsofline (nod[elnod[2]],nod[elnod[3]],x,y) &&
//  	 isnodonlhsofline (nod[elnod[3]],nod[elnod[0]],x,y)     )
//        return (1);
//      break;
//    }
//    case 622:{
//      if ( isnodonlhsof3pcurve (nod[elnod[0]],nod[elnod[3]],nod[elnod[1]],x,y) &&
//  	 isnodonlhsof3pcurve (nod[elnod[1]],nod[elnod[4]],nod[elnod[2]],x,y) &&
//  	 isnodonlhsof3pcurve (nod[elnod[2]],nod[elnod[5]],nod[elnod[0]],x,y)     )
//        return (1);
//      break;
//    }
//    case 822:{
//      if ( isnodonlhsof3pcurve (nod[elnod[0]],nod[elnod[4]],nod[elnod[1]],x,y) &&
//  	 isnodonlhsof3pcurve (nod[elnod[1]],nod[elnod[5]],nod[elnod[2]],x,y) &&
//  	 isnodonlhsof3pcurve (nod[elnod[2]],nod[elnod[6]],nod[elnod[3]],x,y) &&
//  	 isnodonlhsof3pcurve (nod[elnod[3]],nod[elnod[7]],nod[elnod[0]],x,y)     )
//        return (1);
//      break;
//    }
//    //case 413:{ break; }
//    default:{
//      fprintf (stderr,"\n\n wrong dimdegnne in function ispointinel (%s, line %d)",__FILE__,__LINE__);
//      break;
//    }}
//    
//    return (0);
//  }
//  
//  /**
//     Function finds out, whether point (xx,yy) lays on the left hand side of 'line'.
//     Line is defined by two nodes n1,n2.
//     
//     created  3.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  long isnodonlhsofline (gnode &n1,gnode &n2,double xx,double yy)
//  {
//    if (n1.x != n2.x)
//      if ( (n2.x>n1.x ? 1.0:-1.0)*((xx-n1.x)*(n1.y-n2.y)/(n1.x-n2.x)+(n1.y-yy)) <= 1.0e-13 )
//        return (1);
//      else
//        return (0);
//    else
//      if ( (n1.y-n2.y)*(xx-n1.x) > -1.0e-13 )
//        return (1);
//      else
//        return (0);
//  }
//  
//  /**
//     Function finds out, whether point (xx,yy) lays on the left hand side of 'curve'.
//     Curve is defined by three nodes n1,n2,n3.
//     
//     created  3.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  long isnodonlhsof3pcurve (gnode &n1,gnode &n2,gnode &n3,double xx,double yy)
//  {
//    double c;
//  
//    c = (n1.x*(n3.y-n2.y)+n2.x*(n1.y-n3.y)+n3.x*(n2.y-n1.y))/(fabs(n1.x-n3.x)+fabs(n1.y-n3.y));
//    if (-1.0e-13<c && c<1.0e-13)
//      return (isnodonlhsofline (n1,n3,xx,yy));
//    else{
//      // docasne
//      //fprintf (stdout,"\n\n\n hrana je krivocara \n\n\n");
//      // docasne
//      double s,x[3],y[2];
//      
//      c = sqrt((n3.x-n1.x)*(n3.x-n1.x)+(n3.y-n1.y)*(n3.y-n1.y));
//      s = (n3.y-n1.y)/c;
//      c = (n3.x-n1.x)/c;
//      
//      x[0] =  c*(  xx-n1.x) + s*(  yy-n1.y);  y[0] = -s*(  xx-n1.x) + c*(  yy-n1.y);
//      x[1] =  c*(n2.x-n1.x) + s*(n2.y-n1.y);  y[1] = -s*(n2.x-n1.x) + c*(n2.y-n1.y);
//      x[2] =  c*(n3.x-n1.x) + s*(n3.y-n1.y);
//    
//      if ( y[1]*x[0]*x[0]/x[1]/(x[1]-x[2]) + y[1]*x[2]*x[0]/x[1]/(x[2]-x[1]) -y[0] <= 1.0e-13 )
//        return (1);
//      
//      return (0);
//    }
//  }
//  
//  /**
//     Function transforms values from nodes of old mesh into int. points of new mesh.
//     
//     @param gt_old - gtopology of old mesh
//     @param mm_new - mechmat of new mesh
//     @param ipcoord - array of coordinates of new mesh int. points, size is mm_new->tnip x dim
//     @param parentel - array of parent elements of new mesh int. points, size is mt_new->tnip
//     @param nipcomp - number of members (of arrays Mm->ip[].strain , Mm->ip[].stress and Mm->ip[].eqother), which are approximated;
//                      nipcomp[0] <0;Mm->ip[].ncompstr>, nipcomp[1] <0;Mm->ip[].ncompstr>, nipcomp[2] <0;Mm->ip[].ncompother>
//  		    if nipcomp = -1 => all members
//     @param nodvalue - 1D array of nodal values of old mesh, size is gt_old->nn * nvals
//                       value position in array is (val[0] on nod[0] ... val[0] on nod[gt->nn] ... val[nval] on nod[0] ... val[nval] on nod[gt->nn])
//     
//     created  7.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  void transfvalues_ip_indirect (gtopology *gt_old,mechmat *mm_new,double **ipcoord,long *parentel,long nipcomp[],double *nodvalue)
//  {
//    long i,j,maxnchild,nvals;
//    long *nchilds,**childip;
//    double *x,*y,**pointval1,**pointval2,**pointval3;
//    
//    
//    if (nipcomp[0]==-1 || nipcomp[0] > Mm->ip[0].ncompstr  ) nipcomp[0] = Mm->ip[0].ncompstr;
//    if (nipcomp[1]==-1 || nipcomp[1] > Mm->ip[0].ncompstr  ) nipcomp[1] = Mm->ip[0].ncompstr;
//    if (nipcomp[2]==-1 || nipcomp[2] > Mm->ip[0].ncompother) nipcomp[2] = Mm->ip[0].ncompother;
//    
//    nvals = nipcomp[0] + nipcomp[1] + nipcomp[2];
//    
//    nchilds = new long [gt_old->ne];
//    memset (nchilds,0,gt_old->ne*sizeof(long));
//    
//    for (i=0;i<mm_new->tnip;i++)
//      nchilds[parentel[i]]++;
//    
//    maxnchild = 0;
//    childip = new long* [gt_old->ne];
//    for (i=0;i<gt_old->ne;i++){
//      childip[i] = new long [nchilds[i]];
//      if (maxnchild<nchilds[i]) maxnchild = nchilds[i];
//      nchilds[i] = 0;
//    }
//    
//    for (i=0;i<mm_new->tnip;i++)
//      childip[parentel[i]][nchilds[parentel[i]]++] = i;
//    
//    x = new double [maxnchild];
//    y = new double [maxnchild];
//    pointval1 = new double* [maxnchild];
//    pointval2 = new double* [maxnchild];
//    pointval3 = new double* [maxnchild];
//    
//    for (i=0;i<gt_old->ne;i++){
//      for (j=0;j<nchilds[i];j++){
//        x[j] = ipcoord[childip[i][j]][0];
//        y[j] = ipcoord[childip[i][j]][1];
//        
//        if (nipcomp[0]) pointval1[j] = mm_new->ip[childip[i][j]].strain;
//        if (nipcomp[1]) pointval2[j] = mm_new->ip[childip[i][j]].stress;
//        if (nipcomp[2]) pointval3[j] = mm_new->ip[childip[i][j]].eqother;
//      }
//      
//      if (nipcomp[0]) give_valuesinpoints (gt_old,i,nchilds[i],x,y,nipcomp[0],nodvalue                                     ,pointval1,'2');
//      if (nipcomp[1]) give_valuesinpoints (gt_old,i,nchilds[i],x,y,nipcomp[1],nodvalue + gt_old->nn*(nipcomp[0]           ),pointval2,'2');
//      if (nipcomp[2]) give_valuesinpoints (gt_old,i,nchilds[i],x,y,nipcomp[2],nodvalue + gt_old->nn*(nipcomp[0]+nipcomp[1]),pointval3,'2');
//    }
//    
//    
//    delete [] nchilds;
//    for (i=0;i<gt_old->ne;i++)  delete [] childip[i];
//    delete [] childip;
//    delete [] x;
//    delete [] y;
//    for (i=0;i<maxnchild;i++)  pointval1[i] = NULL;  delete [] pointval1;
//    for (i=0;i<maxnchild;i++)  pointval2[i] = NULL;  delete [] pointval2;
//    for (i=0;i<maxnchild;i++)  pointval3[i] = NULL;  delete [] pointval3;
//  }
//  
//  /**
//     Function transforms values from int. points of old mesh into int. points of new mesh.
//     
//     @param gt_old - gtopology of old mesh
//     @param mm_new - mechmat of new mesh
//     @param apcoord - array of coordinates of new mesh int. points, size is mm_new->tnip x dim
//     @param parentel - array of parent elements of new mesh int. points, size is mt_new->tnip
//     @param nipcomp - number of members (of arrays Mm->ip[].strain , Mm->ip[].stress and Mm->ip[].eqother), which are approximated;
//                      nipcomp[0] <0;Mm->ip[].ncompstr>, nipcomp[1] <0;Mm->ip[].ncompstr>, nipcomp[2] <0;Mm->ip[].ncompother>
//  		    if nipcomp = -1 => all members
//     @param mm_old - mechmat of old mesh
//     @param dim - dimension of solved problem
//     
//     created  7.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  void transfvalues_ip_direct (gtopology *gt_old,mechmat *mm_new,const double **apcoord,const long *parentel,long nipcomp[],const mechmat *mm_old,long dim)
//  {
//    long i,j,nvals;
//    double **spvalue,**apvalue;
//    
//    if (nipcomp[0]==-1 || nipcomp[0] > Mm->ip[0].ncompstr  ) nipcomp[0] = Mm->ip[0].ncompstr;
//    if (nipcomp[1]==-1 || nipcomp[1] > Mm->ip[0].ncompstr  ) nipcomp[1] = Mm->ip[0].ncompstr;
//    if (nipcomp[2]==-1 || nipcomp[2] > Mm->ip[0].ncompother) nipcomp[2] = Mm->ip[0].ncompother;
//    
//    nvals = nipcomp[0] + nipcomp[1] + nipcomp[2];
//    
//    spvalue = new double* [gt_old->ne];
//    
//    for (i=0;i<gt_old->ne;i++){
//      switch (gt_old->gelements[i].auxinf){
//      case 312:{
//        spvalue[i] = new double [1*nvals];
//        
//        for (j=0;j<nipcomp[0];j++)  spvalue[i][j                      ] = mm_old->ip[i].strain[j];
//        for (j=0;j<nipcomp[1];j++)  spvalue[i][j+nipcomp[0]           ] = mm_old->ip[i].stress[j];
//        for (j=0;j<nipcomp[2];j++)  spvalue[i][j+nipcomp[0]+nipcomp[1]] = mm_old->ip[i].eqother[j];
//        
//        break;
//      }
//      default:{
//  	fprintf (stderr,"\n\n wrong dimdegnne in function transfvalues_ip_direct (%s, line %d)",__FILE__,__LINE__);
//  	break;
//      }}
//    }
//    
//    //jede jen pro 312 !!!!!!!!!!!!!!!!!!!!
//    apvalue = new double* [mm_new->tnip];
//    for (i=0;i<mm_new->tnip;i++)
//      apvalue[i] = mm_new->ip[i].eqother;    // !!! zde to nereflektuje styl nstra x nstre x nother
//    
//    least_square spr(dim,nvals,gt_old,((Mp->adaptflag & 16) ? 16 : 0));
//    spr.L2_sp2sp (Out,spvalue,mm_new->tnip,parentel,apcoord,apvalue);
//    
//    
//    for (i=0;i<mm_new->tnip;i++) apvalue[i] = NULL;
//    delete [] apvalue;
//  }
//  
//  /**
//     Function transforms displacement from nodes of old mesh into nodes of new mesh by direct interpolation.
//     
//     @param po - problem of old mesh
//     @param pn - problem of new mesh
//     @param lcid - id of load case
//     @param ndofn - number of DOFs in node
//     @param parentel - array of parent (old)elements of new mesh nodes, size is p_new->mt->nn
//     
//     created  7.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  //void transfvalues_nod (problem *po,problem *pn,long lcid,long ndofn,long *parentel)
//  
//  
//  /**
//     Function transforms values from nodes of old mesh into nodes of new mesh.
//     
//     @param po - problem of old mesh
//     @param pn - problem of new mesh
//     @param lcid - id of load case
//     @param dim - dimension of solved problem
//     @param ndofn - number of DOFs in node
//     @param parentel - array of parent (old)elements of new mesh nodes, size is p_new->mt->nn
//     
//     created  7.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  void transfvalues_nod (problem *po,problem *pn,long lcid,long dim,long ndofn,long *parentel,char spr)
//  {
//    long i,j,k,nid,maxnchild;
//    long *nchilds,**childnod;
//    double **spcoord, *r_old, **r_new;
//    vector x,y;
//    
//    // array of deformations of old mesh
//    r_old = new double [ndofn * po->gt->nn];
//    give_rfull (lcid,po,r_old);
//    
//    // allocation
//    nchilds = new long [po->gt->ne];
//    memset (nchilds,0,po->gt->ne*sizeof(long));
//    
//    for (i=0;i<pn->gt->nn;i++)
//      nchilds[parentel[i]]++;
//    
//    childnod = new long* [po->gt->ne];
//    r_new = new double* [po->gt->ne];
//    for (i=0;i<po->gt->ne;i++){
//      childnod[i] = new long [nchilds[i]];
//      r_new[i] = new double [nchilds[i]*ndofn];
//    }
//    
//    if (spr=='y'){
//      spcoord = new double* [po->gt->ne];
//      for (i=0;i<po->gt->ne;i++)  spcoord[i] = new double [nchilds[i]*dim];
//    }
//    else if (spr=='n'){
//      maxnchild = 0;
//      for (i=0;i<po->gt->ne;i++)  if (maxnchild<nchilds[i]) maxnchild = nchilds[i];
//      reallocv (maxnchild,x);  reallocv (maxnchild,y);
//    }
//    else fprintf (stderr,"\n\n unknown spr flag is required in function transfvalues_nod (%s, line %d).\n",__FILE__,__LINE__);
//    
//    // transformation from r_old to r_new
//    memset (nchilds,0,po->gt->ne*sizeof(long));
//    for (i=0;i<pn->gt->nn;i++)
//      childnod[parentel[i]][nchilds[parentel[i]]++] = i;
//    
//    // ********************
//    /*
//    double X,Y,f;
//    
//    for (i=0;i<po->mt->nn;i++){
//      X = po->gt->gnodes[i].x;
//      Y = po->gt->gnodes[i].y;
//      f = -X*X*X*X/1000.0 - X*X*X/10.0 + X*X + 10*X -Y*Y*Y*Y/3.0 - Y*Y*Y/2.0 + 20*Y*Y + 30*Y - 212;
//      
//      r_old[i           ] = f;
//      r_old[i+po->mt->nn] = f;
//    }
//    
//    long hod,min;
//    double sec = clock();
//    */
//    // ********************
//    
//    if (spr=='y'){
//      for (i=0;i<po->gt->ne;i++){
//        for (j=0;j<nchilds[i];j++){
//  	spcoord[i][j*dim  ] = pn->gt->gnodes[childnod[i][j]].x;
//  	spcoord[i][j*dim+1] = pn->gt->gnodes[childnod[i][j]].y;
//  	if (dim==3) spcoord[i][j*dim+2] = pn->gt->gnodes[childnod[i][j]].z;
//        }
//      }
//      
//      least_square spr(dim,ndofn,po->gt,((Mp->adaptflag & 16) ? 16 : 0));
//      spr.L2_nod2sp (Out,nchilds,spcoord,r_old,r_new,'n');
//    }
//    else if (spr=='n'){
//      for (i=0;i<po->gt->ne;i++){
//        for (j=0;j<nchilds[i];j++){
//  	x[j] = pn->gt->gnodes[childnod[i][j]].x;
//  	y[j] = pn->gt->gnodes[childnod[i][j]].y;
//        }      
//        
//        give_valuesinpoints (po->gt,i,nchilds[i],x.a,y.a,ndofn,r_old,&r_new[i],'1');
//      }
//    }
//    else fprintf (stderr,"\n\n unknown spr flag is required in function transfvalues_nod (%s, line %d).\n",__FILE__,__LINE__);
//    
//    // ********************
//    /*
//    sec = (clock() - sec) / (double)CLOCKS_PER_SEC;
//    hod = (long)sec/3600;  sec -= hod*3600;
//    min = (long)sec/60;    sec -= min*60;
//    fprintf (stdout,"\n -----------------------------------");
//    fprintf (stdout,"\n Consumed time by TRANSF %2ld:%02ld:%05.2f",hod,min,sec);
//    fprintf (stdout,"\n -----------------------------------\n");
//    
//    double a,ai,nordr,norr=0.0;
//    vector dr(2*pn->gt->nn),r(pn->gt->nn);
//    
//    for (i=0;i<po->gt->ne;i++)
//      for (j=0;j<nchilds[i];j++){
//        nid = childnod[i][j];
//        X = pn->gt->gnodes[nid].x;
//        Y = pn->gt->gnodes[nid].y;
//        r[nid] = -X*X*X*X/1000.0 - X*X*X/10.0 + X*X + 10*X -Y*Y*Y*Y/3.0 - Y*Y*Y/2.0 + 20*Y*Y + 30*Y - 212;
//        
//        dr[nid+pn->gt->nn] = r[nid] - r_new[i][j*2+1];
//        dr[nid           ] = r[nid] - r_new[i][j*2  ];
//      }
//    
//    sizev (dr,nordr);
//    fprintf (stdout,"\n\n nordr     =  %25.15f \n\n",nordr);
//    fprintf (stdout,"\n\n nordr/nn1 =  %25.15f \n\n",nordr/sqrt(769));
//    fprintf (stdout,"\n\n nordr/nn2 =  %25.15f \n\n",nordr/sqrt(24048));
//    
//    nordr = a = 0.0;
//    for (i=0;i<pn->gt->ne;i++){
//      nordr += Pelt->error(i,dr,ai);
//      a += ai;
//    }
//    
//    nordr = sqrt(nordr/a);
//    fprintf (stdout,"\n\n nordr(A)  =  %25.15f \n\n",nordr);
//    
//    exit (1);
//    
//    
//    for (i=0;i<pn->gt->nn;i++)  norr += r[i];
//    sizev (r,norr); norr *= 2;
//    */
//    // ********************
//    
//    // saving deformations from r_old to p_new->lsrs->lhs
//    for (i=0;i<po->gt->ne;i++)
//      for (j=0;j<nchilds[i];j++){
//        nid = childnod[i][j];
//        for (k=0;k<ndofn;k++)
//  	if (pn->gt->gnodes[nid].cn[k] > 0)
//  	  pn->lsrs->lhs[lcid*pn->lsrs->ndof + pn->gt->gnodes[nid].cn[k] - 1] = r_new[i][j*ndofn+k];
//      }
//    
//    delete [] r_old;
//    delete [] nchilds;
//  }
//  
//  /**
//     Function approximates 'values' from nodes into 'points'.
//     All the point lays in 'element', which is determined by 'gt' and 'eid'.
//     Values in points are returned by array pointval.
//     
//     @param npoints - number of points
//     @param xx - array of x-coordinates of points, size is npoints
//     @param yy - array of y-coordinates of points, size is npoints
//     @param nval - number of values
//     @param nodvalues - 1D array of values in nodes, size is nval*gt->nn, value position in array is
//                        (val[0] on nod[0] ... val[0] on nod[gt->nn] ... val[nval] on nod[0] ... val[nval] on nod[gt->nn])
//     @param pointval - answer = 2D array, size is npoints*nval
//     
//     created  7.12.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  void give_valuesinpoints (gtopology *gt,long eid,long npoints,double *xx,double *yy,long nval,double *nodvalues,double **pointval,char flag)
//  {
//    long i,j,k,nne;
//    double xi,eta,value;
//    vector x,y,nodval;
//    ivector nodes;
//    
//    nne = gt->give_nne (eid);
//    
//    reallocv (nne,x);
//    reallocv (nne,y);
//    reallocv (nne,nodval);
//    reallocv (nne,nodes);
//    
//    gt->give_nodes (eid,nodes);
//    gt->give_node_coord2d (x,y,eid);
//    
//    for (i=0;i<npoints;i++){
//      
//      switch (gt->gelements[eid].auxinf){
//      case 312:{ nc_lin_3_2d          (xx[i],yy[i],x.a,y.a,xi,eta); break; }
//      case 622:{ nc_quad_3_2d (0.00000001,xx[i],yy[i],x.a,y.a,xi,eta); break; }
//      case 412:{ nc_lin_4_2d  (0.00001,xx[i],yy[i],x.a,y.a,xi,eta); break; }
//      case 822:{ nc_quad_4_2d (0.00001,xx[i],yy[i],x.a,y.a,xi,eta); break; }
//      default:{
//        fprintf (stderr,"\n\n unknown nnedegdim is required in function give_valuesinpoints (%s, line %d).\n",__FILE__,__LINE__);
//      }}
//      
//      for (j=0;j<nval;j++){
//        for (k=0;k<nne;k++)
//  	nodval[k] = nodvalues[nodes[k]+j*gt->nn];
//        
//        switch (gt->gelements[eid].auxinf){
//        case 312:{ value = Pelt->approx_nat (xi,eta,nodval);  break; }
//        case 622:{ value = Peqt->approx     (xi,eta,nodval);  break; }
//        case 412:{ value = Pelq->approx     (xi,eta,nodval);  break; }
//        case 822:{ value = Peqq->approx     (xi,eta,nodval);  break; }
//        default:{
//  	fprintf (stderr,"\n\n unknown nnedegdim is required in function give_valuesinpoints (%s, line %d).\n",__FILE__,__LINE__);
//        }}
//        
//        if (flag=='2')
//  	pointval[i][j] = value;
//        else if (flag=='1')
//  	pointval[0][i*nval+j] = value;
//        else ;
//      }
//    }
//  }
//  
/**
   function finds out node(from gt), which is the closest by required point
   
   @param x,y,z - coordinates of point
   
   created  30.1.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
long adjacnode (gtopology *gt,double x,double y,double z)
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
////////////////////       /* termitovo */       ////////////////////////////////////
