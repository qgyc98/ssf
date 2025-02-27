#ifndef MESHTRANSFER_H
#define MESHTRANSFER_H

class problem;
class gtopology;
class mechtop;
class mechmat;
class gnode;
class lhsrhs;

 /* termitovo */
void newmeshread (const char *filename,long lcid);

long give_nodvals_ip (long dim,long nipcomp[],double *&nodvalues);
void give_rfull (long lcid,double *rfull);
void allipcoord (double **ipcoord);

void findout_parentel_ip (gtopology *gt_old,gtopology *gt_new,mechmat *mm_new,mechtop *mt_new,double **ipcoord,long *parentel_ip);
void locintpoints (gtopology *gt_old,gtopology *gt_new,mechtop *mt_new,double **ipcoord,long *passedel1,long *passedel2,long nelheap,long *elheap,long *newelheap,long *susel,long *newsusel,long *parentel);
void findout_parentel_nod (gtopology *gt_old,gtopology *gt_new,long *parentel_nod,long dim);
void locnodes (gtopology *gt_old,gtopology *gt_new,long *parentel,long dim);

long whereispoint (gtopology *gt,long nsusel,long *susel,long *newsusel,long *passedel,double x,double y,char flag);
long ispointinel (gnode *nod,long *elnod,double x,double y,long ndd);
long isnodonlhsofline (gnode &n1,gnode &n2,double xx,double yy);
long isnodonlhsof3pcurve (gnode &n1,gnode &n2,gnode &n3,double xx,double yy);

void transfvalues_ip_indirect (gtopology *gt_old,mechmat *mm_new,double **ipcoord,long *parentel,long nipcomp[],double *nodvalue);
void transfvalues_ip_direct (gtopology *gt_old,mechmat *mm_new,const double **apcoord,const long *parentel,long nipcomp[],const mechmat *mm_old,long dim);

void transfvalues_nod (problem *po,problem *pn,long lcid,long dim,long ndofn,long *parentel,char spr);

void give_valuesinpoints (gtopology *gt,long eid,long npoints,double *xx,double *yy,long nval,double *nodvalues,double **pointval,char flag);

long adjacnode (gtopology *gt,double x,double y,double z);
 /* termitovo */

#endif
