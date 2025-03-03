/*
  File:             homogmat.cpp
  Author:           Jaroslav Kruis, 23.9.2010
  Purpose:          Calculates properties of homogenized material
*/ 

#include "homogmat.h"
#include "globalt.h"
#include "math.h"

homogmat::homogmat (void)
{
  allocm (6,6,dd);
  allocm (3,3,cc);

  hom_mattype = 0;
  hom_mattype_number = 0;
}
homogmat::~homogmat (void)
{
  destrm(dd);
  destrm(cc);
}

/**
   function reads parameters
   
   @param in - input file

   4/8/2017 by TKr
*/
void homogmat::read (XFILE *in)
{
  xfscanf (in,"%ld %ld",&hom_mattype,&hom_mattype_number);
  hom_mattype_number = hom_mattype_number-1;
}


/**
   function prints parameters
   
   @param out - outut file

   4/8/2017 by TKr
*/
void homogmat::print (FILE *out)
{
  fprintf (out,"  %ld %ld",hom_mattype,hom_mattype_number+1);
}



/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   JK, 23.9.2010
*/
void homogmat::matcond (matrix &d,long ri,long ci,long ipp)
{
  long n;
  n = d.n;

  switch (n){
  case 1:{
    matcond1d (d,ri,ci,ipp);//1D
    break;
  }
  case 2:{
    matcond2d (d,ri,ci,ipp);//2D
    break;
  }
  case 3:{
    matcond3d (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    print_err("unknown number of components of conductivity tensor is required",__FILE__,__LINE__,__func__);
  }
  }
}


/**
   function creates conductivity %matrix of the material for 1D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   JK, 23.9.2010
*/
void homogmat::matcond1d (matrix &d,long ri,long ci,long /*ipp*/)
{
  d[0][0] = dd[ri][ci];
}

/**
   function creates conductivity %matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   JK, 23.9.2010
*/
void homogmat::matcond2d (matrix &d,long ri,long ci,long /*ipp*/)
{
  d[0][0] = dd[2*ri+0][2*ci+0];   d[0][1] = dd[2*ri+0][2*ci+1];
  d[1][0] = dd[2*ri+1][2*ci+0];   d[1][1] = dd[2*ri+1][2*ci+1];

  //fprintf(Outt,"\n ipp %ld ri %ld ci %ld  k  %20.15le\n",ipp,ri,ci,d[0][0]);
}

/**
   function creates conductivity %matrix of the material for 3D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   JK, 23.9.2010
*/
void homogmat::matcond3d (matrix &d,long ri,long ci,long /*ipp*/)
{
  d[0][0] = dd[3*ri+0][3*ci+0];   d[0][1] = dd[3*ri+0][3*ci+1];   d[0][2] = dd[3*ri+0][3*ci+2];
  d[1][0] = dd[3*ri+1][3*ci+0];   d[1][1] = dd[3*ri+1][3*ci+1];   d[1][2] = dd[3*ri+1][3*ci+2];
  d[2][0] = dd[3*ri+2][3*ci+0];   d[2][1] = dd[3*ri+2][3*ci+1];   d[2][2] = dd[3*ri+2][3*ci+2];
}


/**
   function creates capacity %matrix of the material

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   JK, 23.9.2010
*/
void homogmat::matcap (double &c,long ri,long ci,long /*ipp*/)
{
  c = cc[ri][ci];

  //fprintf(Outt,"\n ipp %ld ri %ld ci %ld  c  %20.15le\n",ipp,ri,ci,c);
}

/**
   function assembles the conductivity %matrix after homogenization
   
   @param d - array of %matrix entries
   
   JK, 23.9.2010
   completed TKr 12/01/2015
*/
void homogmat::assemble_matrices (double *d,long ntm,long dim)
{
  if (dim==2 && ntm==1){
    dd[0][0]=d[0];   dd[0][1]=d[1];
    dd[1][0]=d[2];   dd[1][1]=d[3];
    
    cc[0][0]=d[4];
  }

  if (dim==3 && ntm==1){

    dd[0][0]=d[0];   dd[0][1]=d[1];   dd[0][2]=d[2];
    dd[1][0]=d[3];   dd[1][1]=d[4];   dd[1][2]=d[5];
    dd[2][0]=d[6];   dd[2][1]=d[7];   dd[2][2]=d[8];
    
    cc[0][0]=d[9];

    /* 
       fprintf (Outt,"\n////////n");
       fprintf (Outt,"\n\n d[0] =  %e",d[0]);
       fprintf (Outt,"\n\n d[1] =  %e",d[1]);
       fprintf (Outt,"\n\n d[2] =  %e",d[2]);
       fprintf (Outt,"\n\n d[3] =  %e",d[3]);
       fprintf (Outt,"\n\n d[4] =  %e",d[4]);
       fprintf (Outt,"\n\n d[5] =  %e",d[5]);
       fprintf (Outt,"\n\n d[6] =  %e",d[6]);
       fprintf (Outt,"\n\n d[7] =  %e",d[7]);
       fprintf (Outt,"\n\n d[8] =  %e",d[8]);
       fprintf (Outt,"\n\n d[9] =  %e",d[9]);
       fprintf (Outt,"\n");
       fprintf (Outt,"\n\n c[0][0] =  %e",cc[0][0]);
    */
    
  }

  if (dim==2 && ntm==2){

    /* fprintf (Outt,"\n////////n");
       fprintf (Outt,"\n\n d[0] =  %e",d[0]);
       fprintf (Outt,"\n\n d[1] =  %e",d[1]);
       fprintf (Outt,"\n\n d[2] =  %e",d[2]);
       fprintf (Outt,"\n\n d[3] =  %e",d[3]);
       fprintf (Outt,"\n\n d[4] =  %e",d[4]);
       fprintf (Outt,"\n\n d[5] =  %e",d[5]);
       fprintf (Outt,"\n\n d[6] =  %e",d[6]);
       fprintf (Outt,"\n\n d[7] =  %e",d[7]);
       fprintf (Outt,"\n\n d[8] =  %e",d[8]);
       fprintf (Outt,"\n\n d[9] =  %e",d[9]);
       fprintf (Outt,"\n\n d[10] =  %e",d[10]);
       fprintf (Outt,"\n\n d[11] =  %e",d[11]);
       fprintf (Outt,"\n\n d[12] =  %e",d[12]);
       fprintf (Outt,"\n\n d[13] =  %e",d[13]);
       fprintf (Outt,"\n\n d[14] =  %e",d[14]);
       fprintf (Outt,"\n\n d[15] =  %e",d[15]);
       fprintf (Outt,"\n\n d[16] =  %e",d[16]);
       fprintf (Outt,"\n\n d[17] =  %e",d[17]);
       fprintf (Outt,"\n\n d[18] =  %e",d[18]);
       fprintf (Outt,"\n\n d[19] =  %e",d[19]);
    */
 
    dd[0][0]=d[0];   dd[0][1]=d[1];   dd[0][2]=d[2];   dd[0][3]=d[3];
    dd[1][0]=d[4];   dd[1][1]=d[5];   dd[1][2]=d[6];   dd[1][3]=d[7];
    dd[2][0]=d[8];   dd[2][1]=d[9];   dd[2][2]=d[10];  dd[2][3]=d[11];
    dd[3][0]=d[12];  dd[3][1]=d[13];  dd[3][2]=d[14];  dd[3][3]=d[15];
    
    cc[0][0]=d[16];    cc[0][1]=d[17];
    cc[1][0]=d[18];    cc[1][1]=d[19];
  }

  if (dim==3 && ntm==2){
    dd[0][0]=d[0];   dd[0][1]=d[1];   dd[0][2]=d[2];   dd[0][3]=d[3];   dd[0][4]=d[4];   dd[0][5]=d[5];
    dd[1][0]=d[6];   dd[1][1]=d[7];   dd[1][2]=d[8];   dd[1][3]=d[9];   dd[1][4]=d[10];  dd[1][5]=d[11];
    dd[2][0]=d[12];  dd[2][1]=d[13];  dd[2][2]=d[14];  dd[2][3]=d[15];  dd[2][4]=d[16];  dd[2][5]=d[17];
    dd[3][0]=d[18];  dd[3][1]=d[19];  dd[3][2]=d[20];  dd[3][3]=d[21];  dd[3][4]=d[22];  dd[3][5]=d[23];
    dd[4][0]=d[24];  dd[4][1]=d[25];  dd[4][2]=d[26];  dd[4][3]=d[27];  dd[4][4]=d[28];  dd[4][5]=d[29];
    dd[5][0]=d[30];  dd[5][1]=d[31];  dd[5][2]=d[32];  dd[5][3]=d[33];  dd[5][4]=d[34];  dd[5][5]=d[35];
    
    cc[0][0]=d[36];  cc[0][1]=d[37];
    cc[1][0]=d[38];  cc[1][1]=d[39];
  }

  /*
  fprintf (stdout,"\n\n homogmat.cpp\n");
  fprintf (stdout,"\n %le %le %le %le",d[0],d[1],d[2],d[3]);
  fprintf (stdout,"\n %le %le %le %le",d[4],d[5],d[6],d[7]);
  fprintf (stdout,"\n %le %le %le %le",d[8],d[9],d[10],d[11]);
  fprintf (stdout,"\n %le %le %le %le",d[12],d[13],d[14],d[15]);

  fprintf (stdout,"\n %le %le %le %le",d[16],d[17],d[18],d[19]);
*/
}


/**
   function computes new transmission coefficient (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param trc     - prescribed transmission coefficient on the boundary
   @param ri      - row index
   @param ci      - column index
   @param nn      - number of node
   @param bc      - type of boundary condition
   @param ipp     - number of first integration point on element

   TKr, 26/09/2011
   corrected by TKr 12/01/2015
*/
double homogmat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  long k;
  double new_trc,w,t;
  new_trc = 0.0;
  
  if(Tp->ntm == 1){
    k=Gtt->give_dof(nn,0);
    if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
    if (k==0)  {t = 0.0;}
    if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
    
    w = 0.0;
    if((ri == 0) && (ci == 0))
      new_trc = get_transmission_transcoeff_tt(w,t,bc,ipp);
  }
  if(Tp->ntm == 2){
    k=Gtt->give_dof(nn,0);
    if (k>0)   {w = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
    if (k==0)  {w = 0.0;}
    if (k<0)   {w = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
    k=Gtt->give_dof(nn,1);
    if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
    if (k==0)  {t = 0.0;}
    if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
    
    if((ri == 0) && (ci == 0))
      new_trc = get_transmission_transcoeff_ww(w,t,bc,ipp);
    if((ri == 0) && (ci == 1))
      new_trc = 0.0;
    
    if((ri == 1) && (ci == 0))
      new_trc = get_transmission_transcoeff_tw(w,t,bc,ipp);
    if((ri == 1) && (ci == 1))
      new_trc = get_transmission_transcoeff_tt(w,t,bc,ipp);
  }

  new_trc = new_trc*trc;

  return (new_trc);
}

/**
   function computes new nodal value (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element

   TKr, 26/09/2011
*/
double homogmat::transmission_nodval(double nodval,double /*trc2*/,long ri,long ci,long nn,long bc,long ipp)
{
  long k;
  double new_nodval,w,t;
  new_nodval = 0.0;

  if(Tp->ntm == 1){
    k=Gtt->give_dof(nn,0);
    if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
    if (k==0)  {t = 0.0;}
    if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
    
    w = 0.0;
    if((ri == 0) && (ci == 0))
      new_nodval = get_transmission_nodval_tt(nodval,w,t,bc,ipp);
  }
  if(Tp->ntm == 2){
    k=Gtt->give_dof(nn,0);
    if (k>0)   {w = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
    if (k==0)  {w = 0.0;}
    if (k<0)   {w = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
    k=Gtt->give_dof(nn,1);
    if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
    if (k==0)  {t = 0.0;}
    if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
    
    if((ri == 0) && (ci == 0))
      new_nodval = get_transmission_nodval_ww(nodval,w,t,bc,ipp);
    if((ri == 0) && (ci == 1))
      new_nodval = 0.0;
    
    if((ri == 1) && (ci == 0))
      new_nodval = get_transmission_nodval_tw(nodval,w,t,bc,ipp);
    if((ri == 1) && (ci == 1))
      new_nodval = get_transmission_nodval_tt(nodval,w,t,bc,ipp);
  }

  return (new_nodval);
}


/**
   function computes flux (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element

   TKr, 26/09/2011
*/
double homogmat::transmission_flux(double nodval,double /*trc2*/,long ri,long ci,long nn,long bc,long ipp)
{
  long k;
  double flux,w,t;
  flux = 0.0;
  
  if(Tp->ntm==1){
    k=Gtt->give_dof(nn,0);
    if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
    if (k==0)  {t = 0.0;}
    if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
    
    w = 0.0;
    if((ri == 0) && (ci == 0))
      flux = get_transmission_flux_tt(nodval,w,t,bc,ipp);
  }
  if(Tp->ntm==2){
    k=Gtt->give_dof(nn,0);
    if (k>0)   {w = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
    if (k==0)  {w = 0.0;}
    if (k<0)   {w = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
    k=Gtt->give_dof(nn,1);
    if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
    if (k==0)  {t = 0.0;}
    if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
    
    if((ri == 0) && (ci == 0))
      flux = get_transmission_flux_ww(nodval,w,t,bc,ipp);
    if((ri == 0) && (ci == 1))
      flux = 0.0;
    
    if((ri == 1) && (ci == 0))
      flux = get_transmission_flux_tw(nodval,w,t,bc,ipp);
    if((ri == 1) && (ci == 1))
      flux = get_transmission_flux_tt(nodval,w,t,bc,ipp);
  }

  return (flux);
}


/**
   function creates transfer coefficient on the boundary for prescribed condition (water vapour pressure)

   @param w   - water content
   @param t   - temperature
   @param bc  - type of boundary condition

   TKr, 26/09/2011
*/
double homogmat::get_transmission_transcoeff_ww(double /*w*/,double /*t*/,long bc,long /*ipp*/)
{
  double trc;
  
  switch (bc){//type of prescribed variable
  case 5:{//water vapour transmission from climatic conditions type 2
    trc = 1.0;
    break;
  }
  case 30:{// relative humidity transmission - it is physically incorrect b.c.; or water vapour transmission
    trc = 1.0;
    break;
  }
  default:{
    print_err("\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  
  return(trc);
}

/**
   function computes the right prescribed value on the boundary for prescribed condition (water vapour pressure)

   @param bv - prescribed value near the boundary
   @param w - actual water conctent on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition

   TKr, 26/09/2011
*/
double homogmat::get_transmission_nodval_ww(double bv,double /*w*/,double /*t*/,long bc,long /*ipp*/)
{  
  double nodval;

  switch (bc){//type of prescribed variable
  case 5:{// water vapour pressure transmission from climatic conditions type 2
    nodval = bv;
    break;
  }
  case 30:{// relative humidity transmission - it is physically incorrect b.c.; or water vapour pressure transmission
    nodval = bv;
    break;
  }
  default:{
    print_err("\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  return(nodval);
}


/**
   function creates flux on the boundary (convective mass transfer)

   @param bv - prescribed value near the boundary
   @param w - water content
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition

   TKr, 26/09/2011
*/
double homogmat::get_transmission_flux_ww(double bv,double w,double /*t*/,long bc,long /*ipp*/)
{
  double flux;

  switch (bc){//type of prescribed variable
  case 5:{// water vapour pressure flux from transmission from climatic conditions type 2
    flux = (bv - w);
    break;
  }
  case 30:{// flux from relative humidity transmission - it is physically incorrect b.c.;or water vapour pressure flux
    flux = (bv - w);
    break;
  }
  default:{
    print_err("\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  return(flux);
}


/**
   function creates transfer coefficient on the boundary for prescribed condition (temperature)

   @param w   - water content
   @param t   - temperature
   @param bc  - type of boundary condition

   TKr, 26/09/2011
*/
double homogmat::get_transmission_transcoeff_tt(double /*w*/,double /*t*/,long bc,long /*ipp*/)
{
  double trc;

  switch (bc){//type of prescribed variable
  case 5:{//heat transmission from climatic conditions type 2
    trc = 1.0;
    break;
  }
  case 30:{//heat transmission
    trc = 1.0;
    break;
  }
  default:{
  case 31:{//heat transmission for testing (and for boundary flux)    
    trc = 0.0;
    break;
  }
    print_err("\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  
  return(trc);
}

/**
   function computes the right prescribed value on the boundary for prescribed condition (temperature)

   @param bv - prescribed value near the boundary
   @param w - actual water conctent on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition

   TKr, 26/09/2011
*/
double homogmat::get_transmission_nodval_tt(double bv,double /*w*/,double t,long bc,long /*ipp*/)
{  
  double nodval;

  switch (bc){//type of prescribed variable
  case 5:{//heat transmission from climatic conditions type 2
    nodval = bv;
    break;
  }
  case 30:{//heat transmission
    nodval = bv;
    break;
  }
  case 31:{//heat transmission for testing (and for boundary flux)
    nodval = (bv - t);
    break;
  }
  default:{
    print_err("\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  return(nodval);
}


/**
   function creates heat flux on the boundary

   @param bv - prescribed value near the boundary
   @param w - water content
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition

   TKr, 26/09/2011
*/
double homogmat::get_transmission_flux_tt(double bv,double /*w*/,double t,long bc,long /*ipp*/)
{
  double flux;

  switch (bc){//type of prescribed variable
  case 5:{//heat flux from heat transmission from climatic conditions type 2
    flux = (bv - t);
    break;
  }
  case 30:{//heat flux from heat transmission
    flux = (bv - t);
    break;
  }
  case 31:{//heat transmission for testing (and for boundary flux)
    flux = (bv - t);
    break;
  }
  default:{
    print_err("\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  return(flux);
}



/**
   function creates transfer coefficient on the boundary for prescribed condition (temperature-moisture)

   @param w   - water content
   @param t   - temperature
   @param bc  - type of boundary condition

   TKr, 26/09/2011
*/
double homogmat::get_transmission_transcoeff_tw(double /*w*/,double /*t*/,long bc,long /*ipp*/)
{
  double trc;

  switch (bc){//type of prescribed variable
  case 5:{//heat transmission from water vapour transmission from climatic conditions type 2
    trc = 1.0;
    break;
  }
  case 30:{//heat transmission from relative humidity transmission - physically incorrect; or from water vapour transmission
    trc = 1.0;
    break;
  }
  case 31:{//heat transmission for testing (and for boundary flux)    
    trc = 0.0;
    break;
  }
  default:{
    print_err("\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  
  return(trc);
}

/**
   function computes the right prescribed value on the boundary for prescribed condition (temperature-moisture)

   @param bv - prescribed value near the boundary
   @param w - actual water conctent on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition

   TKr, 26/09/2011
*/
double homogmat::get_transmission_nodval_tw(double bv,double /*w*/,double /*t*/,long bc,long /*ipp*/)
{  
  double nodval;

  switch (bc){//type of prescribed variable
  case 5:{//heat transmission from water vapour transmission from climatic conditions type 2
    nodval = bv;
    break;
  }
  case 30:{//heat transmission from relative humidity transmission - physically incorrect; or from water vapour transmission
    nodval = bv;
    break;
  }
  case 31:{//heat transmission for testing (and for boundary flux)
    nodval = 0.0;
    break;
  }
  default:{
    print_err("\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  return(nodval);
}


/**
   function creates heat flux on the boundary from moisture

   @param bv - prescribed value near the boundary
   @param w - water content
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition

   TKr, 26/09/2011
*/
double homogmat::get_transmission_flux_tw(double bv,double w,double /*t*/,long bc,long /*ipp*/)
{
  double flux;

  switch (bc){//type of prescribed variable
  case 5:{//heat flux from water vapour transmission from climatic conditions type 2
    flux = (bv - w);
    break;
  }
  case 30:{//heat flux from water vapour transmission
    flux = (bv - w);
    break;
  }
  case 31:{//heat transmission for testing (and for boundary flux)
    flux = 0.0;
    break;
  }
  default:{
    print_err("\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  return(flux);
}



/**
   function computes all variables in nodes
   @param compother - number of other components
   @param ipp - first integration point on element
   @param pc - capillary pressure on actual node
   @param pg - gas pressure on actual node
   @param t - temperature on actual node
*/

double homogmat::get_othervalue(long compother,double rh,double t,long /*ipp*/)
{
  double other,tk;
  
  if(Tp->ntm == 1){
    switch (compother){
    case 0:{//temperature
      other = t;
      break;
    }
    default:{
      print_err("\n\n This component type is not implemented\n",__FILE__,__LINE__,__func__);
    }
    }
  }
  if(Tp->ntm == 2){
    switch (compother){
    case 0:{//relative humidity
      if (hom_mattype == 155)
	other = rh;
      if (hom_mattype == 180){
	if(Tm->moisth == NULL){
	  print_err("Basic material for homogenization is not defined in macro problem",__FILE__,__LINE__,__func__);
	  abort();
	}
	tk = Tm->moisth[hom_mattype_number].teptr.getval(0);
	
	//fprintf(Outt,"Tm->moisth[hom_mattype_number].teptr.f = %e\n",Tm->moisth[hom_mattype_number].teptr.f);
	//fprintf(Outt,"hom_mattype_number = %ld\n",hom_mattype_number);
	//fprintf(Outt,"hom_mattype = %ld\n",hom_mattype);
	//fprintf(Outt,"tk = %e\n",tk);
	
	other = rh/(exp(23.5771 - 4042.9/(tk - 37.58)));
      }
      break;
    }
    case 1:{//temperature
      other = t - 273.15;
      break;
    }
    case 2:{//temporarilly
      other = 0.0;
      break;
    }
    case 3:{//temporarilly
      other = 0.0;
      break;
    }
    case 4:{//temporarilly
      other = 0.0;
      break;
    }
    case 5:{//temporarilly
      other = 0.0;
      break;
    }
    case 6:{//temporarilly
      other = 0.0;
      break;
    }
    case 7:{//temporarilly
      other = 0.0;
      break;
    }
    case 8:{//temporarilly
      other = 0.0;
      break;
    }
    case 9:{//temporarilly
      other = 0.0;
      break;
    }
    case 10:{//temporarilly
      other = 0.0;
      break;
    }
    default:{
      print_err("\n\n This component type is not implemented\n",__FILE__,__LINE__,__func__);
    }
    }
  }

  return (other);
}

/**
   function prints names of all variables in nodes
   @param out - output file
   @param compother - number of other components
   
   TKr, 26/09/2011
*/
void homogmat::print_othervalue_name(FILE *out,long compother)
{
  if(Tp->ntm == 1){
    switch (compother){
    case 0:{//Temperature
      fprintf (out,"Temperature (K)             ");
      break;
    }
    default:{
      print_err("\n\n This component type is not implemented\n",__FILE__,__LINE__,__func__);
    }
    }
  }
  if(Tp->ntm == 2){
    switch (compother){
    case 0:{//Relative humidity
      fprintf (out,"Relative humidity ()       ");
      break;
    }
    case 1:{//Temperature
      fprintf (out,"Temperature (K)             ");
      break;
    }
    default:{
      print_err("\n\n This component type is not implemented\n",__FILE__,__LINE__,__func__);
    }
    }
  }
}

/**
   The function returns ordered dof names of primary unknowns 
   required by the model.
   
   @param dofname   - array of uknown name for particular nodal dofs (output)
                      dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
   @param ntm       - number of transported media = number of nodal dof = length of array dofname
   
   JK, 23. 4. 2014
*/
void homogmat::give_dof_names(namevart *dofname, long ntm)
{
  if (ntm < 1)
  {
    print_err("the model defines %ld unknowns but number of transported media is %ld", 
              __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }

  if(ntm == 1)
    dofname[0] = trf_temperature;

  if(ntm == 2){
    //  these names are provisional
    if (hom_mattype == 180){
      dofname[0] = trf_press_water_vapor;
      dofname[1] = trf_temperature;
    }
    else{
      dofname[0] = trf_rel_humidity;
      dofname[1] = trf_temperature;      
    }
  }
}

/**
   The funtion marks required non-transport quantities in the array antq.
   
   @param antq - array with flags for used material types
                 antq[i] = 1 => quantity type nontransquant(i+1) is required
                 antq[i] = 0 => quantity type nontransquant(i+1) is not required

   @return The function does not return anything, but it may change content of antq array.
   
   TKr 01/08/2022 according to TKo
*/
void homogmat::give_reqntq(long */*antq*/)
{
}
