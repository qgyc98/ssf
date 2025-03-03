#include <mpi.h>
#include "fixnodesel.h"

/**
   Function round off the number in agreement with mathematical conventions
   Function return double
   @param number - input number for round off
   
   JB
 **/

double round_it(double number)
{
  double fractpart,intpart;
  fractpart = modf (number,&intpart);
  if(fractpart >= 0.5){
    number=ceil(number);
  }
  else{
    number=floor(number);
  }
  return(number);
}


/**
   Function check if the input integer is prime number or not
   @param number - input number for prime number detection
   Function return logical value - true/false
   true - input number is prime number
   false - input number is not prime number
   
   JB
 **/
bool isPrime(int number) 
{
  long i,max;
 
  if (number < 2) return false;
  if (number == 2) return true;
  if (number % 2 == 0) return false;
  
  max = (long) sqrt(number);
  
  for (i = 3; i <= max; i += 2) {
    if (number % i == 0) {
      return false;
    }
  }
  return true;
}

/**
   Function compute number of combination
   @param n 
   @param k
   ( n )
   ( k )
   Output is @double format - number of combination

   JB
 **/


long compute_number_of_combination(long n, long k)
{
  long i;
  long denom,nom,diff;
  long combination;
  
  diff = n - k;
  nom = 1;
  for(i = diff+1; i < n+1; i++){
    nom *= i;
  }
  denom = 1;
  for(i = 2; i < k+1; i++){
    denom *= i;
  }
  combination = nom/denom;
  
  return(combination);
}


/**
   Constructor of the @class fixnodesel
   @param np - number of processors
   @param mr - myrank
   @param nd - number of subdomain
   @param meshd - mesh description
   @param nameproc - name of processor
   @param nl - length of name of processor
   @param long mes - message printing
   JB
 **/

fixnodesel::fixnodesel(int np,int mr,long nd,meshdescription meshd,long* domproces,char *nameproc, int nl,long mes)
{
  mespr=mes;

  //  number of processors
  nproc=np;
  //  my rank
  myrank=mr;
  //  number of subdomain
  ndom=nd;
  //domain processor corespomdence
  domproc = domproces;
  
  // name of processor
  strcpy (procName,nameproc);
  nameLength=nl;
  
  //  mesh description
  md = meshd;
  
  // array with coordinates of fictitious subdomain
  fictdom = NULL;
  // array with interface nodeidentification on subdomain
  nodeidentif = NULL;
  // array with interface nodeidentification on master
  coarseidentif = NULL;
  // array with numbers of boundary nodes adjacent to boundary node
  coarsenadjac = NULL;
  // list of boundary nodes adjacent to boundary node
  coarseadjac = NULL;
  // array with identification of nodes on boundary curves
  curidentif=NULL;
  // array with identification of nodes on boundary surfaces
  surfidentif=NULL;
  // identifiation if boundary node lie on the surface
  surfnod=NULL;
  // list of boundary nodes adjacent to boundary node on bounbdary curve
  curadjac=NULL;
  // array with numbers of boundary nodes adjacent to boundary node on bounbdary curve
  curnadjac=NULL;
  // identification of boundary node lie on the curve
  curnod=NULL;
  // array with start nodes of boundary curves
  start=NULL;
  // array with end nodes of boundary curves
  end=NULL;
  // list of members of boundary curves
  members=NULL;
  // number of edges of boundary curves
  nedges=NULL;
  // number of members of boundary curves
  nmembers=NULL;
  // list of position of fixing nodes defined by user
  userdefnod=NULL;
  // number of boundary surfaces
  nsurf = 0;
  // list of members of boundary surfaces
  surfmembers = NULL;
  // numbers of members of boundary surfaces
  nsurfmembers = NULL;
  // 
  automember = NULL;
  // global subdomain - boundary node correspondence
  glinkdom = NULL;
  // global node - local node correspondence
  glinknod = NULL;
  // local subdomain - node correspondence
  loclinkdom = NULL;
  // array with boundary node indentification on master processor
  coarseidentif = NULL;
  
  // array with centre of gravity of boundary surfaces
  surfcenters = NULL;
  // number of centre of gravity of boundary surfaces
  nsurfcentres = NULL;  
  // array of indicators if the boundary node lie on the boundary surfaces
  surfnodpoint = NULL;
  // array with marks of boundary surface nodes 
  surfnodmark = NULL;
  // number of boundary surfaces    
  nsurf = 0;
  // number of user specified position of fixing nodes
  nuserdefnod = -1;
  // spatial dimension of problem
  dim = -1;
  // selection strategy for minimal number algorithm 
  minSelStrat = -1;
  // type of condesation of fixing nodes
  condfixing = nomethod;
  // method of condensation of boudnary curves
  typecondcur = notype;
  // number of fixing nodes added onto boundary curves
  nmembercur = 0;
  // method of condensation of boudnary surfaces
  typecondsurf = notype;
  // number of fixing nodes added onto boundary surfaces
  nmembersurf = 0;
  // place of condensation of fixing nodes
  methodcondcor = 0;
  // size of fictitious subdomain
  sizefictdom = -1.0;
  // array with numbers of boundary nodes adjacent to boundary node
  nadjacboundnod = NULL;
  // list of boundary nodes adjacent to boundary node
  adjacboundnod = NULL;

  // tolerance for comparsion of two values of double type
  tol = 1e-10;
  // ration between side of fixtitious subdomain
  ratio = NULL;
  
}

/**
   Destructor of @class fixnodesel
 
**/

fixnodesel::~fixnodesel()
{
  long i,j;


  if(fictdom != NULL){
    delete []fictdom;
  }
  if( nodeidentif != NULL){
    delete []nodeidentif;
  }
  if( userdefnod != NULL){
    delete []userdefnod;
  }
  
  if(myrank == 0){
    if(coarseadjac != NULL){
      for(i = 0; i < tnbn; i++){
	delete []coarseadjac[i];
      }
    }
    if(coarsenadjac != NULL){
      delete []coarsenadjac;  
    }
    
    if(coarseidentif != NULL){
      delete []coarseidentif;
    }
    if(surfidentif != NULL){
      delete []surfidentif;
    }
     if(curidentif != NULL){
       delete []curidentif;
     }
     if(curadjac != NULL){
       for(j = 0; j < ncurnod; j++){
	 delete []curadjac[i];
       }
       delete []curadjac;
     }
     if(curnadjac != NULL){
       delete []curnadjac;
     }
     if(curnod != NULL){
       delete []curnod;
     }
    if( start != NULL){
      delete []start;
    }
    if( end != NULL){
      delete []end;
    }
    if( members != NULL){
      for(j = 0; j < ncurves; j++){
	delete []members[j];
      }
      delete []members;
    }
    if( nmembers != NULL){
      delete []nmembers; 
    }
    if(nedges != NULL){
      delete []nedges;
    }
    
    if(surfnod != NULL){
      delete []surfnod;
    }
    if(surfnodpoint != NULL){
      // delete []surfnodpoint;
    }
  }   
  
  if(surfmembers != NULL){
    for(i = 0; i < nsurf; i++){
      delete surfmembers[i];
    }
    delete []surfmembers;
  }
  if(nsurfmembers != NULL){
    delete []nsurfmembers;
  }
  
  if(glinkdom != NULL){
    for(i = 0; i < tnbn; i++){
      delete []glinkdom[i];
    }
    delete []glinkdom;
  }
  if(glinknod != NULL){
    for(i = 0; i < tnbn; i++){
      delete []glinknod[i];
    }
    delete []glinknod;
  }
  
  if(loclinkdom != NULL){
    for(i = 0; i < nbn; i++){
      delete []loclinkdom[i];
    }
    delete []loclinkdom;
  }
    
  if(automember != NULL){
    delete []automember;
  }
  if(surfcenters != NULL){
    for(i = 0; i < nsurf; i++){
      delete []surfcenters[i];
    }
    delete []surfcenters;
  }
  if(nsurfcentres != NULL){
    delete []nsurfcentres;
  }  
  
  
  if(surfnodmark != NULL){
    for(i = 0; i < nsurf; i++){
      delete []surfnodmark[i];
    }
    delete []surfnodmark;
  }
 
  delete []ratio;
  if(mespr == 1) fprintf(out,"\n\n\nFIXNODESEL CLASS - end\n");  
}

/** 
    Function initialize @class fixnodesel
    @param topol - pointer to gtopology class
    @param part - pointer to partop class
    @param outfile - pointer to outputfile
          
    JB
 **/


void fixnodesel::initiate(gtopology *topol,partop *part,FILE *outfile)
{
  // out file
  out = outfile;
  // general topology class
  top = topol;
  // partop class
  ptop = part;
  //  number of nodes on subdomain
  nn = top->nn;
  //  number of elements on subdomain
  ne = top->ne;
  
  if (myrank==0){
    //nbndod - list with numbers of boundary nodes on subdomains on master
    nbnd = ptop->nbnd;
    //list of multiplicity of boudary nodes on master
    multip = ptop->bmultip;
    //array with coarse numbers of boundary nodes on master
    // cnbn[i][j] = k the j-th node on the i-th subdomain has coarse number k
    cnbn = ptop->icnbnmas;
    //total number of bouundary nodes on master
    tnbn = ptop->tnbn;
  }
  //list of nodal multiplocity of nodes on subdomain
  nodmultip = ptop->nodmultip;
  //number of boundary nodes on subdomain
  nbn = ptop->nbn;
  // list of local numbers of boundary nodes on subdomain
  lnbn = ptop->lnbndom;
  //list with coarse numbers of nodes on subdomain
  lgnbn = ptop->icnbn;
  //maximum number of boudnary nodes on subdomain
  maxnbn = ptop->maxnbn;
  //number of internal nodes on subdomain
  nin = ptop->nin;
  //list of local numbers of internal nodes on subdomain
  lnin= ptop->lnindom;
  
  if(mespr == 1) fprintf(out,"\n\n\nFIXNODESEL CLASS - start\n");  

}

/**
   Function checks array ltg. If there are fixing nodes then class fixnodesel
   is finished. If there are not fixnig nodes class fixnodesel selects appropriate
   fixing nodes. Fuction is called from class psolver - psolver.cpp.
   @param ltg - array with global - local correspondence
   @param nnd - number of nodes on subdomain
   @param out - pointer to output file for logs
   Function return bool type - false if fixing nodes are found in ltg array
                             - true  if fixing nodes are not found in ltg array
   JB
 */
bool fixnodesel::check_ltg(long *ltg,long nnd,FILE *out)
{
  long neg,pos,one,i;
  switch (md){
  case all_nodes:{
    return false;
    break;
  }
  case bound_nodes:{
    neg = 0;
    pos = 0;
    one = 0;
    for (i = 0; i < nnd; i++ ){
      if (ltg[i] <= -2){
	neg++;
      }
      if (ltg[i] > -1){
	pos++;
      }
      if (ltg[i] == -1){
	one++;
      }
    }
    //fprintf(out,"Cheking neg=%ld pos=%ld one=%ld\n",neg,pos,one);
    //if ((neg+pos+one) == nn){
    //fprintf(out,"Cheking OK\n");
    //}
    // fixing nodes are found
    if ( neg > 0){
      if(mespr == 1) fprintf(out,"\n\n\nFixing nodes are found in ltg array\n");
      return false;
    }
    // fixing nodes must be found
    else{
      if(mespr == 1) fprintf(out,"\n\n\nFixing nodes must be found by function fixing_detection\n");
      return true;
    }
    break;
  }
  case neg_bound_nodes:{
    return false;
    break;
  }
  default:{
    
  }
  }
  fflush(out);
  return false; // added due to compiler warning
}


/**
   function recognizes spatial dimension of the mesh
   @output dim - spatial dimension
   29.5.2007, JB
*/

void fixnodesel::give_whole_dim()
{
  long *elem;
  long i;
 
  elem = new long[ne];
  
  for(i = 0; i < ne; i++){
    elem[i] =  top->give_whole_dim(i);
  }
  for(i = 1; i < ne; i++){
    if(elem[i] != elem[i-1]){
      par_print_err(myrank,"different spatial dimension is found in the mesh,function fixing_detection can not be used\n", __FILE__, __LINE__, __func__);
    }
    else{
      dim = elem[i];
    }
  }
  
  delete []elem;
}




/**
   Function computes lenght of vector which is determined by nodes a and b
   @param a - start node
   @param b - end node
   u=(u1,u2,u3)
   u1 = xb-xa
   u2 = yb-ya
   u3 = zb-za
   @output u - length of vector

   JB
**/

double fixnodesel::compute_length_of_vector(long a,long b)
{
  double u1,u2,u3,x1,x2,y1,y2,z1,z2,u;
    
  //node 1 
  x1 = top->gnodes[a].x;
  y1 = top->gnodes[a].y;
  z1 = top->gnodes[a].z;
  // node 2
  x2 = top->gnodes[b].x;
  y2 = top->gnodes[b].y;
  z2 = top->gnodes[b].z;
  // vecotor u - nodes 1 and 2
  u1 = x2 - x1;
  //fprintf(out,"u1 = %le\n",u1);
  u2 = y2 - y1;
  //fprintf(out,"u2 = %le\n",u2);
  u3 = z2 - z1;
  //  magnitude of vectoru u
  u = u1*u1+u2*u2+u3*u3;
  u=sqrt(u);
  return(u);
}

/**
   Function computes angle between two vectors  which are defined by nodes a, b and c
   vector u = b - a
   vector v = c - a
   Output is variable double angle.
   
   @param top - pointer to the sequential general topology
   @param a - the first node of vector u and vector v
   @param b - the second node of vector u
   @param c - the second node of vector v
   @param out - output file (used for auxiliary output)
   
   JB 24.09.2009
**/

double fixnodesel::compute_angle_of_vector_a_and_b(long a,long b,long c)
{
  double u1,u2,u3,v1,v2,v3,x1,x2,x3,y1,y2,y3,z1,z2,z3,scal,u,v,cosalpha,angle;
  
  //node 1 
  x1 = top->gnodes[a].x;
  y1 = top->gnodes[a].y;
  z1 = top->gnodes[a].z;
  // node 2
  x2 = top->gnodes[b].x;
  y2 = top->gnodes[b].y;
  z2 = top->gnodes[b].z;
  // vecotor u - nodes 1 and 2
  u1 = x2 - x1;
  //fprintf(out,"u1 = %le\n",u1);
  u2 = y2 - y1;
  //fprintf(out,"u2 = %le\n",u2);
  u3 = z2 - z1;
  //fprintf(out,"u3 = %le\n",u3);
  // node 3
  x3 = top->gnodes[c].x;
  y3 = top->gnodes[c].y;
  z3 = top->gnodes[c].z; 
  // vector v - nodes 1 and 3
  v1 = x3 - x1;
  //fprintf(out,"v1 = %le\n",v1);
  v2 = y3 - y1;
  //fprintf(out,"v2 = %le\n",v2);
  v3 = z3 - z1;
  //fprintf(out,"v3 = %le\n",v3);
  // scalar product
  scal = u1*v1+u2*v2+u3*v3;
  //fprintf(out,"scal = %le\n",scal);
  //  magnitude of vector  u
  u = u1*u1+u2*u2+u3*u3;
  //fprintf(out,"u = %le\n",u);
  u = sqrt(u);
  //fprintf(out,"u = %le\n",u);
  // magnitude of vector  v
  v = v1*v1+v2*v2+v3*v3;
  //fprintf(out,"v = %le\n",v);
  v = sqrt(v);
  //fprintf(out,"v = %le\n",v);
  // angle
  cosalpha = scal/(u*v);
  
  if(cosalpha > 1.0){
    cosalpha = 1.0;
  }
  if(cosalpha < -1.0){
    cosalpha = -1.0;
  }
  angle = acos(cosalpha);
  return(angle);

}

/**
   Function compute area of triangle which is defined by nodes a, b and c.
   Computation is based on cross product of vectors
   Output is variable double area.
   
   @param top - pointer to the sequential general topology
   @param a - the first node of triangle
   @param b - the second node of triangle
   @param c - the third node of triangle
   @param out - output file (used for auxiliary output)

   
   JB 24.09.2009
**/

double fixnodesel::compute_area_of_triangle(long a,long b,long c)
{
  //fprintf(out,"compute area of triangle\n");
  double x1,x2,x3,y1,y2,y3,z1,z2,z3,u1,u2,u3,v1,v2,v3,w1,w2,w3,w,area;
  //double pr,dr;
  //double v1,v2,alpha,sinalpha;
  //   double alpha,sinalpha,area;
  //   pr=compute_length_of_vector(a,b);
  //   dr=compute_length_of_vector(a,c);
  //   alpha=compute_angle_of_vector_a_and_b(a,b,c);
  //   alpha=acos(alpha);
  //   sinalpha=sin(alpha);
  //   area=(pr*dr*sinalpha)/2;
  

  //node 1 
  x1 = top->gnodes[a].x;
  y1 = top->gnodes[a].y;
  z1 = top->gnodes[a].z;
  // node 2
  x2 = top->gnodes[b].x;
  y2 = top->gnodes[b].y;
  z2 = top->gnodes[b].z;
  // vecotor u - nodes 1 and 2
  u1 = x2 - x1;
  //fprintf(out,"u1 = %le\n",u1);
  u2 = y2 - y1;
  //frintf(out,"u2 = %le\n",u2);
  u3 = z2 - z1;
  //fprintf(out,"u3 = %le\n",u3);
  // node 3
  x3 = top->gnodes[c].x;
  y3 = top->gnodes[c].y;
  z3 = top->gnodes[c].z; 
  // vector v - nodes 1 and 3
  v1 = x3 - x1;
  //fprintf(out,"v1 = %le\n",v1);
  v2 = y3 - y1;
  //fprintf(out,"v2 = %le\n",v2);
  v3 = z3 - z1;
  //fprintf(out,"v3 = %le\n",v3);
  
  //cross product c=uxv
  w1=u2*v3-u3*v2;
  w2=u3*v1-u1*v3;
  w3=u1*v2-u2*v1;
  // magnitude of vector  c
  w = w1*w1+w2*w2+w3*w3;
  //fprintf(out,"w = %le\n",w);
  area=w/2;
  
  
//   area = x1*y2*z3+y1*z2*x3+z1*x2*y3-z1*y2*x3-y1*x2*z3-x1*z2*y3;
  
  return(area);
}


/** 
    Function check geometrical condition of selection of fixing nodes
    JB
**/
void fixnodesel::check_triangle()
{
  long i,j,k;
  double *length,*angle;
  long v;
  long *vcomb,*ver,*pointver;
  long stop;

  //fprintf(out,"check triangle\n");
  
  vcomb = new long[3];
  
  v = 0 ;
  for(i = 0; i < nbn; i++){
    if(nodeidentif[i] == 3){
      v++;
    }
  }
  ver = new long[v];
  pointver = new long[v];
  j = 0 ;
  for(i = 0; i < nbn; i++){
    if(nodeidentif[i] == 3){
      ver[j] = lnbn[i];
      pointver[j] = i;
      j++;
    }
  }
  
  length = new double[3];
  angle = new double[3];
  for(i = 0; i < v; i++){
    stop = 0;
    if(nodeidentif[pointver[i]] == 3){
      vcomb[0] = ver[i];
      for(j = i+1; j < v; j++){
	stop = 0;
	if(nodeidentif[pointver[j]] == 3){
	  vcomb[1] = ver[j];
	  length[2]= compute_length_of_vector(vcomb[0],vcomb[1]);
	  for(k = j+1; k < v; k++){
	    if(nodeidentif[pointver[k]] == 3){
	      vcomb[2] = ver[k];
	      //fprintf(out,"uzly : %ld %ld %ld\n",vcomb[0],vcomb[1],vcomb[2]);
	      length[1] = compute_length_of_vector(vcomb[0],vcomb[2]);
	      length[0] = compute_length_of_vector(vcomb[1],vcomb[2]);
	      angle[0] = compute_angle_of_vector_a_and_b(vcomb[0],vcomb[1],vcomb[2]);
	      length[0] = compute_length_of_vector(vcomb[1],vcomb[2]);
	      angle[1] = compute_angle_of_vector_a_and_b(vcomb[1],vcomb[0],vcomb[2]);
	      angle[2] = compute_angle_of_vector_a_and_b(vcomb[2],vcomb[0],vcomb[1]);
	      //test uhlu
	      //if(angle[0]+angle[1]+angle[2] > 3.1416 || angle[0]+angle[1]+angle[2] < 3.1414){
		//fprintf(out,"divne uhly\n");
	      //}
	      //test print
	      //fprintf(out,"delky : %le %le %le\n",length[0],length[1],length[2]);
	      //fprintf(out,"uhly : %lf %lf %lf\n",angle[0]*180/3.141592653589793238462643,angle[1]*180/3.141592653589793238462643,angle[2]*180/3.141592653589793238462643);
	      // vyhodnoceni
	      // kontrola primky
	      //if(fabs(cos(angle[0])) > 0.99){
	      if((angle[0] < 0.034906585 || angle[0] > 3.1066861) || (angle[1] < 0.034906585 || angle[1] > 3.1066861) || (angle[2] < 0.034906585 || angle[2] > 3.1066861)){
	      // lezi na primce
		if(length[1] >= length[2]){
		  nodeidentif[pointver[j]]=2;
		  stop = 1;
		  //fprintf(out,"uzly lezi na primce, bude vymazan uzel %ld\n",lnbndom[pointver[j]]);
		}
		else{
		  nodeidentif[pointver[k]]=2;
		  //fprintf(out,"uzly lezi na primce, bude vymazan uzel %ld\n",lnbndom[pointver[k]]);
		}
	      }
	      
	      // toto neni dotazeno do konce
	      // kontrola uhlu 
	      if(angle[0] < (3.141592653589793238462643*5)/180 && angle[1] < (3.141592653589793238462643*5)/180){
		nodeidentif[pointver[k]]=2;
		//fprintf(out,"male uhly 0 a 1 bude vymazan uzel %ld\n",lnbndom[pointver[k]]);
	      }
	      if(angle[0] < (3.141592653589793238462643*5)/180 && angle[2] < (3.141592653589793238462643*5)/180 ){
		nodeidentif[pointver[j]]=2;
		//fprintf(out,"male uhly 0 a 2 bude vymazan uzel %ld\n",lnbndom[pointver[j]]);
		
	      }
	      if(angle[1] < (3.141592653589793238462643*5)/180 && angle[2] < (3.141592653589793238462643*5)/180){
		nodeidentif[pointver[i]]=2;
		//fprintf(out,"male uhly 1 a 2 bude vymazan uzel %ld\n",lnbndom[pointver[i]]);
	      }
	      
	      // kontrola delek
	      if(length[0]< 0.05*sizefictdom ){
		if(length[1] > 0.05*sizefictdom && length[1] > length[2]){
		  nodeidentif[pointver[j]]=2;
		  //fprintf(out,"mala vzdalenost mezi uzly bude vymazan uzel %ld\n",lnbndom[pointver[j]]);
		  stop = 1;
		}
		if(length[2] > 0.05*sizefictdom && length[2] > length[1]){
		  //fprintf(out,"mala vzdalenost mezi uzly bude vymazan uzel %ld\n",lnbndom[pointver[k]]);
		  nodeidentif[pointver[k]]=2;
		}
	      }
	      if(length[1]< 0.05*sizefictdom ){
		if(length[0] > 0.05*sizefictdom && length[0] > length[2]){
		  //fprintf(out,"mala vzdalenost mezi uzly, bude vymazan uzel %ld\n",lnbndom[pointver[i]]);
		  nodeidentif[pointver[i]]=2;
		  stop = 2;
		}
		if(length[2] > 0.05*sizefictdom && length[2] > length[0]){
		  nodeidentif[pointver[k]]=2;
		  //fprintf(out,"mala vzdalenost mezi uzly, bude vymazan uzel %ld\n",lnbndom[pointver[k]]);
		}
	      }
	      if(length[2]< 0.05*sizefictdom ){
		if(length[0] > 0.05*sizefictdom && length[0] > length[1]){
		  nodeidentif[pointver[i]]=2;
		  //fprintf(out,"mala vzdalenost mezi uzly, bude vymazan uzel %ld\n",lnbndom[pointver[i]]);
		  stop = 2;
		}
		if(length[1] > 0.05*sizefictdom && length[1] > length[0]){
		  nodeidentif[pointver[j]]=2;
		  //fprintf(out,"mala vzdalenost mezi uzly, bude vymazan uzel %ld\n",lnbndom[pointver[j]]);
		  stop = 1;
		}
	      }
	      if(stop == 1){
		break;
	      }
	    }
	  }
	  if(stop == 2){
	    break;
	  }
	}
      }
    }
  }

  j = 0;
  if(mespr == 1){
    for(i = 0; i < nbn; i++){
      //fprintf(out,"uzel %ld nodeidentif %ld\n",lnbndom[i],nodeidentif[i]);
      if(nodeidentif[i] == 3){
	fprintf(out,"Node %ld is fixing nodes\n",lnbn[i]);
	j++;
      }
    }
    fprintf(out,"Number of fixing node is %ld, after check_triangle\n",j);
  }
  delete []vcomb;
  delete []ver;
  delete []pointver;
  
  delete []length;
  delete []angle;
}

/*
  function establishes fictitious subdomain around_it real subdomain
  fictitious subdomain is rectangle in two dimension
  fictitious subdomain is cube in three dimension

  fictdom[0] = k - min_x coordinates
  fictdom[1] = k - min_y coordinates
  fictdom[2] = k - min_z coordinates
  fictdom[3] = k - max_x coordinates
  fictdom[4] = k - max_y coordinates
  fictdom[5] = k - max_z coordinates

  10.02.2009 JB
 */

void fixnodesel::set_fictitious_subdomain()
{
  long i;
  double  miny,minx,minz,maxy,maxx,maxz,x,y,z;

  if(fictdom != NULL){
    delete []fictdom;
  }
  fictdom = new double[6]; 
  //miny=DBL_MAX;
  //minx=DBL_MAX;
  //minz=DBL_MAX;
  //maxy=DBL_MIN;
  //maxx=DBL_MIN;
  //maxz=DBL_MIN;
  minx = top->gnodes[0].x;
  miny = top->gnodes[0].y;
  minz = top->gnodes[0].z;
  maxx = top->gnodes[0].x;
  maxy = top->gnodes[0].y;
  maxz = top->gnodes[0].z;
  for(i = 0; i < nn; i++){
    x = top->gnodes[i].x;
    y = top->gnodes[i].y;
    z = top->gnodes[i].z;	
    //fprintf(out,"%ld %lf %lf %lf\n",i,x,y,z);
    // min_x
    if(x <= minx){
      minx = x;
      fictdom[0] = x;
    }
    // min_y
    if(y <= miny){
       miny = y;
      fictdom[1] = y;
    }
    // min_z
    if(z <= minz){
      minz = z;
      fictdom[2] = z;
    }
    // max_x
    if(x >= maxx){
      maxx = x;
      fictdom[3] = x;
    }
    // max_y
    if(y >= maxy){
      maxy = y;
      fictdom[4] = y;
    }
    // max_z
    if(z >= maxz){
      maxz = z;
      fictdom[5] = z;
    }
  }
  if(mespr == 1) fprintf(out,"fictdom: %lf %lf %lf %lf %lf %lf\n",fictdom[0],fictdom[1],fictdom[2],fictdom[3],fictdom[4],fictdom[5]);
  fflush(out);
}


/**
   function computes distance between nodes with minimal and maximal coordinates of fictitious subdomain
   10.02.2009 JB
   lx/ly
   ratio = 1 square
   ratio > 1
   y
    -------------
   |             |
    -------------   x
   ratio < 1 
   y
    --
   |  |
   |  |
   |  |
    --  x
	      
 */

void fixnodesel::compute_size_of_fictitious_subdomain()
{
  double vetcx,vetcy,vetcz,sizevect;
  //size_fictdom
  vetcx = fictdom[3] - fictdom[0];
  vetcy = fictdom[4] - fictdom[1];
  vetcz = fictdom[5] - fictdom[2];
  //fprintf(out,"velikost fictdom %lf %lf %lf\n",vetc_x,vetc_y,vetc_z);
  sizevect = vetcx*vetcx + vetcy*vetcy + vetcz*vetcz;
  sizefictdom = sqrt(sizevect);
  if(mespr == 1) fprintf(out,"size of fictdom %lf\n",sizefictdom);
  lx=fictdom[3]-fictdom[0];
  ly=fictdom[4]-fictdom[1];
  lz=fictdom[5]-fictdom[2];

  
  ratio = new double[3];
  if( ly != 0.0){ 
    ratio[0]=lx/ly;
  }
  else{
    ratio[0] = 0.0;
  }
  if( lz != 0.0){ 
    ratio[1]=lx/lz;
  }
  else{
    ratio[1] = 0.0;
  }
  if( lz != 0.0){
    ratio[2]=ly/ly;
  }
  else{
    ratio[2] = 0.0;
  }
  if(mespr == 1) fprintf(out,"size of fictdom: x side %lf, y side %lf,z side %lf\n",lx,ly,lz);
  if(mespr == 1) fprintf(out,"ratio: xy %lf, xz  %lf,yz %lf\n",ratio[0],ratio[1],ratio[2]);
  fflush(out);
}

/**
   Function compute statistics for selection of minimal number of fixing nodes
   @output is minSelStrat
   minSelStrat = 1 - choice is based on statistics
   minSelStrat = 2 - choice is based on graph theory
   JB
 **/

void fixnodesel::compute_statistics_of_multiplicity()
{
  long i,j,k,l;
  long buffsize;
  long *buff;
  MPI_Status stat;
  
  if(mespr == 1) fprintf(out,"\n\n\nStatistics:\n");
  
  maxmultip = 0;
  for(i = 0; i < nn; i++){
    if(nodmultip[i] > maxmultip){
      maxmultip = nodmultip[i];
    }
  }
  
  if(myrank == 0){
    //  array containing number of nodes on subdomains
    nstatdom = new long [nproc];
    statdom = new long*[nproc];
    //  master contribution
    j=domproc[0];
    nstatdom[j]=maxmultip;
    statdom[j] = new long[nstatdom[j]];
    for (i = 1; i < nproc; i++){
      MPI_Recv (&k,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //  slave contributions
      j=domproc[stat.MPI_TAG];
      if (maxmultip < k){
	maxmultip = k;
      }
      nstatdom[j] = k;
      statdom[j] = new long[maxmultip];
    }
    
    if(mespr == 1){
      fprintf(out,"\n\n\nMaximal multiplicity  %ld\n",maxmultip);
      for (i = 0; i < nproc; i++){
	fprintf(out,"\n\nDomain %ld has maximal multiplicity %ld\n",i,nstatdom[i]);
      }
    }
    
    
    for (i=1;i<nproc;i++){
      MPI_Send (&maxmultip,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Send (&maxmultip,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (&maxmultip,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  buffsize = maxmultip;
  buff = new long[buffsize];
  
  for(j = 0; j < maxmultip; j++){
    buff[j] = 0;
    for(i = 0; i < nn; i++){
      if(nodmultip[i] == j+1){
	buff[j]++;
      }
    }
    //fprintf(out,"Pocet uzlu s multiplicitou %ld je %ld\n",j+1,buff[j]);
  }
  
  if(myrank == 0){
    k=domproc[0];
    for(j = 0; j < maxmultip; j++){
      statdom[k][j] = buff[j];
    }
    
    for (i = 1; i < nproc; i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        
      //  slave contributions
      k=domproc[stat.MPI_TAG];
      for(j = 0; j < maxmultip; j++){
      statdom[k][j] = buff[j];
      }
    }
    
    if(mespr == 1){
      for (i = 0; i < nproc; i++){
	fprintf(out,"\n\nDomain %ld:\n",i);
	for(j= 0; j < maxmultip; j++){
	  fprintf(out,"Number of nodes with multiplicity %ld is %ld\n",j+1,statdom[i][j]);
	}
      }
    }
    
    long estMinCornod;
    l = 0;
    estMinCornod = nproc*3;
    if(mespr == 1) fprintf(out,"\n\nEstimate minimal number of  fixing nodes in whole problem %ld\n",estMinCornod);
    for (i = 0; i < nproc; i++){
      k = 0;
      for(j= 1; j < maxmultip; j++){
	if(statdom[i][j] < estMinCornod){
	  k++;
	}
      }
      if(k == 1){
	l++;
      }
    }
    if(l != 0){
      minSelStrat = 1;
    }
    else{
      minSelStrat = 2;
    }
    if(mespr == 1) fprintf(out,"\n\nselected strategy: %ld\n",minSelStrat);
    for (i = 1; i < nproc; i++){
      MPI_Send (&minSelStrat,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (&minSelStrat,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete []buff;
}

/**
   Function creates global node - domain correspondence and global node - local node corespondence
   @output glinknod
   glinknod has tnbn components
   glinknod[i] has multip[i] components
   glinknod[i][j] = k - 
   @output glinkdom
   glinkdom has tnbn components
   glinkdom[i] has multip[i] components
   glinkdom[i][j] = k - 
   JB

**/
void fixnodesel::create_link_dom_nod()
{
  long i,j,m,k,l,max;
  MPI_Status stat;
  long *buff,buffsize;
  
  
  if(myrank == 0){
    
    // vazby
    glinknod = new long*[tnbn];
    glinkdom = new long*[tnbn];
    
    for(i = 0; i < tnbn; i++){
      glinknod[i] = new long[multip[i]];
      glinkdom[i] = new long[multip[i]];
    }
    
    for(i = 0; i < tnbn; i++){
      m = 0;
      for(j = 0; j < nproc; j++){
	for(k = 0; k < nbnd[j]; k++){
	 if(cnbn[j][k] == i){
	   glinknod[i][m] = k;
	   glinkdom[i][m] = j;
	   m++;
	   break;
	 }
       }
      }
    }
    max = 0;
    for(i = 0; i < nproc; i++){
      m = 0;
      for(j = 0; j < tnbn; j++){
	for(k = 0; k < multip[j]; k++){
	  if(glinkdom[j][k] == i){
	    m += multip[j] - 1;
	  }
	}
      }
      if(m > max){
	max = m;
      }
    }
    //fprintf(out,"max je %ld \n",max);
    
    buffsize = max;
    for(i = 1; i < nproc; i++){
      MPI_Send (&buffsize,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
//     for(i = 0; i < tnbn; i++){
//       fprintf(out,"%ld glinkdom na ",i);
//       for(j = 0; j < multip[i]; j++){
// 	fprintf(out,"  %ld",glinkdom[i][j]);
//       }
//       fprintf(out,"\n");
//     }
//     for(i = 0; i < tnbn; i++){
//       fprintf(out,"%ld glinknod na ",i);
//       for(j = 0; j < multip[i]; j++){
// 	fprintf(out,"  %ld",glinknod[i][j]);
//       }
//       fprintf(out,"\n");
//     }
  }
  // ********
  // slaves
  // ********
  else{
    MPI_Recv (&buffsize,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  max = buffsize;

  buff = new long[buffsize];

  
  if(myrank == 0){
    
    for(i = 1; i < nproc; i++){
      l = 0;
      for(j = 0; j < tnbn; j++){
	for(k =0; k < multip[j]; k++){
	  if(glinkdom[j][k] == i){
	    for(m = 0; m < k; m++){
	      buff[l]=glinkdom[j][m];
	      l++;
	    }
	    for(m = k+1; m < multip[j]; m++){
	      buff[l]=glinkdom[j][m];
	      l++;
	    }
	    break;
	  }
	}
      }
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    l = 0;
    for(j = 0; j < tnbn; j++){
      for(k =0; k < multip[j]; k++){
	if(glinkdom[j][k] == 0){
	  for(m = 0; m < k; m++){
	    buff[l]=glinkdom[j][m];
	      l++;
	  }
	  for(m = k+1; m < multip[j]; m++){
	    buff[l]=glinkdom[j][m];
	    l++;
	  }
	  break;
	}
      }
    }
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  
  loclinkdom = new long*[nbn];
  l = 0;
  for(i = 0; i < nbn; i++){
    loclinkdom[i] = new long[nodmultip[lnbn[i]]-1];
    for(j = 0; j < nodmultip[lnbn[i]]-1; j++){
      loclinkdom[i][j]=buff[l];
      l++;
    }
  }
  delete []buff;
  // for(i = 0; i < nbn; i++){
//     fprintf(out,"node %ld loclinkdom:",i);
//     for(j = 0; j < nodmultip[lnbn[i]]-1; j++){
//       fprintf(out,"   %ld",loclinkdom[i][j]);
//     }
//     fprintf(out,"\n");
//   }
  

}

/**
   Function creates subdomain boundary graph for 2D meshes
   @output nadjacboundnod - array with numbers of nodes adjacent to nodes
   nadjacboundnod has nbn components
   adjacboundnod[i] = j - the i-th boundary node has j adjacent nodes
   @output adjacnodes
   adjacboundnod has nbn components
   adjacboundnod[i] has nodmultip[i] components
   adjacboundnod[i][j] = k the j-th adjacent node to i-th node has number k
   JB
**/
void fixnodesel::create_subdom_graph_2d()
{
  long i,j,k,l,m;
  long ned,nned,auxnedges;
  long *auxnned,*dupledges,*aux,*locbound;
  long **auxedgenodes;
  
  // number of all edges in mesh
  auxnedges = 0;
  // pomocne pole
  aux = new long[3];
  for(i = 0; i < ne; i++){
    // number of edge on i-th element - gtopology
    ned = top->give_ned(i);
    //fprintf(out,"ned %ld\n",ned);
    nned = top->give_nned(i);
    //fprintf(out,"nned %ld\n",nned);
    for(j = 0; j < ned; j++ ){
      top->give_edge_nodes(i,j,aux);
      m = 0;
      for(l = 0; l < nned; l++){
	if(nodmultip[aux[l]] > 1){
	  m++;
	}
      }
      if(m > 1){
	auxnedges++;
      }
    }
  }
  //fprintf(out,"pocet hran %ld\n",auxnedges);
  
  // plneni pole s uzky na hranach
  // plneni s pole s poctem uzlu na hrane
  auxnned = new long[auxnedges];
  auxedgenodes = new long*[auxnedges];
  dupledges = new long[auxnedges];
  // pomocne pole
  //aux = new long[3];
  k = 0;
  for(i = 0; i < ne; i++){
    ned = top->give_ned(i);
    nned = top->give_nned(i);
    for(j = 0; j < ned; j++ ){
      top->give_edge_nodes(i,j,aux);
      m = 0;
      for(l = 0; l < nned; l++){
	if(nodmultip[aux[l]] > 1){
	  m++;
	}
      }
      if(m > 1){
	auxnned[k] = nned;
	auxedgenodes[k] = new long[nned];
	for(l = 0; l < nned; l++){
	  auxedgenodes[k][l] = aux[l];
	}
	// trideni pole s uzly na hrane(od nejmensiho po nejvetsi)
	if(nned == 2){
	  if(auxedgenodes[k][0] > auxedgenodes[k][1]){
	    l = auxedgenodes[k][1];
	    auxedgenodes[k][1] = auxedgenodes[k][0];
	    auxedgenodes[k][0] = l;
	  }
	}
	if(nned == 3){
	  if(auxedgenodes[k][0] > auxedgenodes[k][2]){
	    l = auxedgenodes[k][2];
	    auxedgenodes[k][2] = auxedgenodes[k][0];
	    auxedgenodes[k][0] = l;
	  }
	  if(auxedgenodes[k][0] > auxedgenodes[k][1]){
	    l = auxedgenodes[k][1];
	      auxedgenodes[k][1] = auxedgenodes[k][0];
	      auxedgenodes[k][0] = l;
	    }
	  if(auxedgenodes[k][1] > auxedgenodes[k][2]){
	    l = auxedgenodes[k][2];
	    auxedgenodes[k][2] = auxedgenodes[k][1];
	    auxedgenodes[k][1] = l;
	  }
	}
	dupledges[k] = 0;
	k++;
      }
    }
  }
  delete []aux;
  //auxnedges = k;
  //fprintf(out,"pocet hran %ld\n",auxnedges);
  
  
  // kontrolni tisk
  //for(i = 0; i < auxnedges; i++){
  //fprintf(out,"hrana %ld ma uzly:",i);
  //for(j = 0; j < auxnned[i]; j++ ){
  //  fprintf(out,"   %ld",auxedgenodes[i][j]);
  //}
  //fprintf(out,"\n");
  //}
  
  // hledani duplikovanych hran na hranici
  for(i = 0; i < auxnedges; i++){
    for(j = i+1; j < auxnedges; j++){
      if(auxedgenodes[i][0] == auxedgenodes[j][0]){
	if(auxnned[i] == 2 && auxnned[j] == 2){
	  if(auxedgenodes[i][1] == auxedgenodes[j][1]){
	    dupledges[j]++;
	    dupledges[i]++;
	  }
	}
	if(auxnned[i] == 3 && auxnned[i] == 3){
	  if(auxedgenodes[i][1] == auxedgenodes[j][1]){
	    if(auxedgenodes[i][2] == auxedgenodes[j][2]){
	      dupledges[j]++;
	      dupledges[i]++;
	    }
	  }
	}
      }
    }
  }
  
  // kontrolni tisk
  //for(i = 0; i < auxnedges; i++){
  //fprintf(out,"%ld dupledges %ld\n",i,dupledges[i]);
  //}
  
  locbound = new long[nn];
  m = 0;
  for(i = 0; i < nn; i++){
    if(nodmultip[i] > 1){
      locbound[i] = m;
      m++;
    }
    else{
      locbound[i] = -1;
    }
  }
  
  // kontrolni tisk
  //     for(i = 0; i < nbn; i++){
  //     fprintf(out,"%ld lnbn %ld\n",i,lnbn[i]);
  //     }
  //     for(i = 0; i < nn; i++){
  //     fprintf(out,"%ld locbound %ld\n",i,locbound[i]);
  //     }
  
  // renumbering
  //if(adjacboundnod != NULL){
  //for(j = 0; j < nadjacboundnod[i]; j++){
  //delete []adjacboundnod[i];
  //}
  //delete []adjacboundnod;
  //}
  //if(nadjacboundnod != NULL){
  //delete []nadjacboundnod;
  //}
  
    
  nadjacboundnod = new long[nbn];
  for(i = 0; i < nbn; i++){
    nadjacboundnod[i] = 0;
    for(j = 0; j < auxnedges; j++){
      if(dupledges[j] == 0) {
	for(k = 0; k < auxnned[j]; k++){
	  if(auxedgenodes[j][k] == lnbn[i]){
	    auxedgenodes[j][k] = i;
	    nadjacboundnod[i]++;
	  }
	}
      }
    } 
  }
  
  // kontrolni tisk
  //for(i = 0; i < auxnedges; i++){
  //fprintf(out,"hrana %ld ma uzly:",i);
  ///for(j = 0; j < auxnned[i]; j++ ){
  // fprintf(out,"   %ld",auxedgenodes[i][j]);
  //}
  //fprintf(out,"\n");
  //}
  //for(i = 0; i < nbn; i++){
  //fprintf(out,"%ld\n",nadjacboundnod[i]);
  //}  
  
  // sestaveni sousedu hranicnich uzlu (pouze na hranici) pro pripadne sestaveni krivek
  adjacboundnod = new long*[nbn];
  for(i = 0; i < nbn; i++){
    adjacboundnod[i] = new long[nadjacboundnod[i]];
    m = 0;
    for(j = 0; j < auxnedges; j++){
      if(dupledges[j] == 0) {
	if(auxnned[j] == 2){
	  for(k = 0; k < auxnned[j]; k++){
	    if(auxedgenodes[j][k] == i){
	      if(k == 0){
		adjacboundnod[i][m] = auxedgenodes[j][1];
		m++;
	      }
	      if(k == 1){
		adjacboundnod[i][m] = auxedgenodes[j][0];
		m++;
	      }
	    }
	  }
	}
	if(auxnned[j] == 3){
	  for(k = 0; k < auxnned[j]; k++){
	    if(auxedgenodes[j][k] == lnbn[i]){
	      if(k == 0){
		adjacboundnod[i][m] = auxedgenodes[j][1];
		m++;
	      }
	      if(k == 1){
		adjacboundnod[i][m] = auxedgenodes[j][0];
		m++;
		adjacboundnod[i][m] = auxedgenodes[j][2];
		m++;
	      }
	      if(k == 2){
		adjacboundnod[i][m] = auxedgenodes[j][1];
		m++;
	      }
	    }
	  }
	}
      }
    }
  }
  // kontrolni tisk
  //for(i = 0; i < nbn; i++){
  // fprintf(out,"uzel %ld ma  %ld sousedu: ",i,nadjacboundnod[i]);
  //for(j = 0; j < nadjacboundnod[i]; j++){
  //fprintf(out,"   %ld",adjacboundnod[i][j]);
  //}
  //fprintf(out,"\n");
  //}
  //for(i = 0; i < nbn; i++){
  // fprintf(out,"uzel %ld ma lgnbn %ld\n",i,lgnbn[i]);
  //}
  
  for(i = 0; i < nbn; i++){
    for(j = 0; j < nadjacboundnod[i]; j++){
      k = adjacboundnod[i][j];
      adjacboundnod[i][j] = lgnbn[k];
    }
  }
  
  
  //sorting 
  for(i = 0; i < nbn; i++){
    if(nadjacboundnod[i] == 2){
      if(adjacboundnod[i][0] > adjacboundnod[i][1]){
	j = adjacboundnod[i][1];
	adjacboundnod[i][1] = adjacboundnod[i][0];
	adjacboundnod[i][0] = j;
      }
    }
  }
  
  
  // kontrolni tisk
 //  for(i = 0; i < nbn; i++){
//     fprintf(out,"uzel %ld ma  %ld sousedu: ",i,nadjacboundnod[i]);
//     for(j = 0; j < nadjacboundnod[i]; j++){
//       fprintf(out,"   %ld",adjacboundnod[i][j]);
//     }
//     fprintf(out,"\n");
//   }
  
  // mazani poli
  for(i = 0; i < auxnedges; i++){
    delete []auxedgenodes[i];
  }
  delete []auxnned;
  delete []auxedgenodes;
  delete []dupledges;
  delete []locbound;

}

/**
   Function select minimal number of fixing nodes for 2D meshes
   Selection is based on graph theory and is described in dissertation JB
   @output coarseidentif - array with identification of boundaty nodes
   coarseidentif = 2 - remaining interface node
   coarseidentif = 3 - fixing node
   coarseidentif has tnbn components
   JB
 **/
void fixnodesel::select_minimal_number_2d()
{
  long i,auxnv;
  
  if(myrank == 0){
    // identifikace uzlu na domene
    if (coarseidentif != NULL)
      delete []coarseidentif;
    coarseidentif = new long [tnbn];
    
    auxnv = 0;
    //fprintf(out,"tnbn je %ld\n\n",tnbn);
    for(i = 0; i < tnbn; i++){
      //origin or terminus of master graph
      if(multip[i] == 2 && coarsenadjac[i] == 1){
	coarseidentif[i] = 3;
      }
      if(multip[i] > 2 && coarsenadjac[i] > 2){
	coarseidentif[i] = 3;
      }
      if(multip[i] == 2 && coarsenadjac[i] == 2){
	coarseidentif[i] = 2;
      }
      if(coarseidentif[i] == 3){
	auxnv++;
      }
    }
    if(mespr == 1) fprintf(out,"Number of potencial fixing nodes in whole problem is %ld\n\n",auxnv);
    fflush(out);
    
    
    if(mespr == 1){
      for (i = 0; i < tnbn; i++){
	if(coarseidentif[i] == 3) fprintf(out,"\n node number %ld is fixing nodes node",i);
      }
      fflush(out);
      fprintf(out,"\n\n");
    }
  }
}

/**
   Function check minimal number of selected fixing nodes for 2D meshes
   JB
 **/
void fixnodesel::check_minimal_number_2d()
{
  long a,b,i,j,k,l,m,auxnv;
  long buffsize,*buff;
  MPI_Status stat;
  long *nvdom;
  long **stor;
  
  buffsize = maxnbn+1;
  buff = new long[buffsize];
  if(myrank == 0){
    nvdom = new long[nproc];
    for(i = 0; i < nproc; i++){
      nvdom[i] = 0;
    }
    
    for(i = 0; i < tnbn; i++){
      if(coarseidentif[i] == 3){
	for(j = 0; j < multip[i]; j++){
	  k = glinkdom[i][j];
	  nvdom[k]++;
	}
      }
    }
    // for(i = 0; i < nproc; i++){
    //       fprintf(out,"domena %ld nvdom %ld\n\n",i,nvdom[i]);    
    //     }
    
    //slave
    for(i = 1; i < nproc; i++){
      for(j = 0; j < nbnd[i]; j++){
	k = cnbn[i][j];
	buff[j] = coarseidentif[k];
      }
      buff[maxnbn] = nvdom[i];
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    //master
    for(j = 0; j < nbnd[0]; j++){
      k = cnbn[0][j];
      buff[j] = coarseidentif[k];
    }
    buff[maxnbn] = nvdom[0];
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  nodeidentif = new long[nbn];
  
  auxnv = 0;
  for(i = 0; i < nbn; i++){
    nodeidentif[i] = buff[i];
    if(nodeidentif[i] == 3){
      auxnv++;
      if(mespr == 1) fprintf(out,"local number of fixng node: %ld\n\n",lnbn[i]);
    }
  }
  
  
  if(mespr == 1) fprintf(out,"Number of potential fixing nodes on subdomain is %ld\n\n",buff[maxnbn]);
  fflush(out);
  

  if(buff[maxnbn] < 2){
    auxnv = buff[maxnbn];
    select_fixing_nodes_geom_2d(auxnv);
    buff[maxnbn] = -2;
  }
  else{
    if(condfixing == nocondconer){
      check_triangle();
      auxnv = 0;
      for(i = 0; i < nbn; i++){
	if(nodeidentif[i] == 3){
	  auxnv++;
	}
      }
      if(auxnv != buff[maxnbn]){
	buff[maxnbn] = -1*auxnv;
      }
      if(auxnv < 2){
	select_fixing_nodes_geom_2d(auxnv);
	buff[maxnbn] = -2;
      }
    }
  }
  
  if(buff[maxnbn] < 0){
    for(i = 0; i < nbn; i++){
      buff[i] = nodeidentif[i]; 
    }
  }
  if(mespr == 1) fprintf(out,"Number of potential fixing nodes %ld on subdomain %d\n\n",abs(buff[maxnbn]),myrank);
  
  
  
  if(myrank == 0){
    //master
    stor = new long*[nproc];
    //fprintf(out,"domena 0 nvdom stare %ld nove %ld\n\n",nvdom[0],buff[maxnbn]);
    if (buff[maxnbn] < 0 ){
      nvdom[0] = -1*buff[maxnbn];
      stor[0] = new long [nbnd[0]];
      for(j = 0; j < nbnd[0]; j++){
      	stor[0][j] = buff[j];
      }
    }
    else{
      stor[0]= new long [nbnd[0]];
      for(j = 0; j < nbnd[0]; j++){
      	stor[0][j] = buff[j];
      }
    }
    
    //slave
    for(i = 1; i < nproc; i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);  
      k=domproc[stat.MPI_TAG];
      //fprintf(out,"domena %ld nvdom stare %ld nove %ld\n\n",k,nvdom[k],buff[maxnbn]);
      if (buff[maxnbn] < 0 ){
	nvdom[k] = -1*buff[maxnbn];
	stor[k] = new long [nbnd[k]];
	for(j = 0; j < nbnd[k]; j++){
	  stor[k][j] = buff[j];
	}
      }
      else{
	stor[k] = new long [nbnd[k]];
	for(j = 0; j < nbnd[k]; j++){
	  stor[k][j] = buff[j];
	}
      }
    }
    
    
    for(i = 0; i < tnbn; i++){
      // controlling
      m = 0;
      l = 0;
      if (coarseidentif[i] == 3){
	a = 1;
	b = 0;
      }
      else{
	a = 0;
	b = 1;
      }	
      for(j = 0; j < multip[i]; j++){
	if(coarseidentif[i] != stor[glinkdom[i][j]][glinknod[i][j]]){
	  if(stor[glinkdom[i][j]][glinknod[i][j]] == 3){
	    a++;
	  }
	  if(stor[glinkdom[i][j]][glinknod[i][j]] == 2){
	    b++;
	  }
	  m++;
	}
	if(nvdom[glinkdom[i][j]] <= 2){
	  l++;
	}
      }
      if( m != 0){
	//fprintf(out,"m = %ld a= %ld b = %ld\n",m,a,b);
	if(a >= b){
	  if(mespr == 1) fprintf(out,"prepis %ld z 2 na 3\n",i);
	  coarseidentif[i] = 3;
	  for(j = 0; j < multip[i]; j++){
	    if(stor[glinkdom[i][j]][glinknod[i][j]] == 2){
	      nvdom[glinkdom[i][j]]++;
	    }
	  }
	}
	else{
	  if(l == 0){
	    if(mespr == 1) fprintf(out,"prepis %ld z 3 na 2\n",i);
	    coarseidentif[i] = 2;
	    for(j = 0; j < multip[i]; j++){
	      if(stor[glinkdom[i][j]][glinknod[i][j]] == 3){
		nvdom[glinkdom[i][j]]--;
	      }
	    }
	  }
	  else{
	    
	  }
	}
      }
    }
    if(mespr == 1){
      for(i = 0; i < nproc; i++){
	fprintf(out,"\n\n\n Number of fixing nodes on domain %ld is %ld\n",i,nvdom[i]);
      }
    }
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);  

  delete []buff;
  
  //info
  if(myrank == 0){
    m = 0;
    for(i = 0; i < nproc; i++){
      nvdom[i] = 0;
    }
    for(i = 0; i < tnbn; i++){
      if(coarseidentif[i] == 3){
	m++;
	for(j = 0; j < multip[i]; j++){
	  k = glinkdom[i][j];
	  nvdom[k]++;
	}
      }
    }
    if(mespr == 1){
      fprintf(out,"\n\n\n Total number of fixing nodes in whole problem is %ld\n",m);
      for(i = 0; i < nproc; i++){
	fprintf(out,"\n\n\n Number of fixing nodes on domain %ld is %ld\n",i,nvdom[i]);
      }
    }
    
    delete []nvdom;
    for(i = 0; i < nproc; i++){
      if(stor[i] != NULL){
	delete []stor[i];
      }
    }
    delete []stor;
  }

}

/**
   Function selects fixing nodes based on geometrical condition of 2D mesh
   JB
 **/
void fixnodesel::select_fixing_nodes_geom_2d(long auxnv)
{

  long i,j,max;
  long *ver;
  double maxradius,radius;
  long *vcomb;

  //double t1= clock();
  switch(auxnv){
  case 0:
    ver = new long[2];
   //  max = 0;
    //     for(i = 0; i < nbn; i++){
    //       if(max < lgnbn[i]){
    // 	max = lgnbn[i];
    // 	ver[0] = i;
    //       }
    //     }
    //     //fprintf(out,"max je %ld\n\n",max);
    //     //fprintf(out,"prvni je %ld\n\n",lnbn[ver[0]]);
    
    // inicialization of pseudorandom numbers
    srand((unsigned int) time(NULL));
    j = rand()%(nbn-1);

    ver[0] = j;
    
    vcomb = new long[2];
    vcomb[0] = lnbn[ver[0]];
    maxradius = 0.0;
    for(i = 0; i < nbn; i++){
      vcomb[1] = lnbn[i];
      radius = compute_length_of_vector(vcomb[0],vcomb[1]);
      //fprintf(out,"radius je %lf\n\n",radius);
      if(radius >= maxradius){
	maxradius = radius;
	j = i;
      }
    }
    ver[0] = j;
    // 	fprintf(out,"novy je %ld\n\n",lnbn[ver[0]]);
    // 	fprintf(out,"maxradius je %lf\n\n",maxradius);
    
    vcomb[0] = lnbn[ver[0]];
    maxradius = 0.0;
    for(i = 0; i < nbn; i++){
      vcomb[1] = lnbn[i];
      radius = compute_length_of_vector(vcomb[0],vcomb[1]);
      if(radius >= maxradius){
	j = i;
	maxradius = radius;
      }
    }
    ver[1] = j;
    nodeidentif[ver[0]]=3;
    nodeidentif[ver[1]]=3;
    vcomb[0] = lnbn[ver[0]];
    vcomb[1] = lnbn[ver[1]];
    // 	fprintf(out,"potencialni vrcholy jsou %ld %ld\n\n",lnbn[ver[0]],lnbn[ver[1]]);
    // 	fprintf(out,"maxradius je %lf\n\n",maxradius);
    
    delete []ver;
    delete []vcomb;
    
    // 	for(i = 0; i < nbn; i++){
    // 	  if(nodeidentif[i] == 3){
    // 	    fprintf(out,"%ld je fixing\n\n",lnbn[i]);
    // 	  }
    // 	}
    //double t2=  clock();
    //fprintf(out,"cas hledani je %lf\n\n",(t2-t1)/(double)CLOCKS_PER_SEC);    
    break;
  case 1:
    ver = new long[2];
    max = 0;
    for(i = 0; i < nbn; i++){
      if(nodeidentif[i] == 3){
	ver[0] = i;
      }
    }
    //fprintf(out,"max je %ld\n\n",max);
    //fprintf(out,"prvni je %ld\n\n",lnbndom[ver[0]]);
    
    vcomb = new long[2];
    vcomb[0] = lnbn[ver[0]];
    maxradius = 0.0;
    for(i = 0; i < nbn; i++){
      vcomb[1] = lnbn[i];
      radius = compute_length_of_vector(vcomb[0],vcomb[1]);
      //fprintf(out,"radius je %lf\n\n",radius);
      if(radius >= maxradius){
	maxradius = radius;
	j = i;
      }
    }
    ver[1] = j;
    //fprintf(out,"dalsi je %ld\n\n",lnbndom[ver[1]]);
    //fprintf(out,"maxradius je %lf\n\n",maxradius);
    //fprintf(out,"potencialni vrcholy jsou %ld %ld\n\n",lnbndom[ver[0]],lnbndom[ver[1]]);
    
    nodeidentif[ver[1]]=3;
    vcomb[0] = lnbn[ver[0]];
    vcomb[1] = lnbn[ver[1]];
    
    delete []ver;
    delete []vcomb;
    
    // 	for(i = 0; i < nbn; i++){
    // 	  if(nodeidentif[i] == 3){
    // 	    fprintf(out,"%ld je fixing\n\n",lnbn[i]);
    // 	  }
    // 	}
    break;
  }
  
}

/**
Function creates subdomain boundary graph for 3D meshes
   @output nadjacboundnod - array with numbers of nodes adjacent to nodes
   nadjacboundnod has nbn components
   adjacboundnod[i] = j - the i-th boundary node has j adjacent nodes
   @output adjacnodes
   adjacboundnod has nbn components
   adjacboundnod[i] has nodmultip[i] components
   adjacboundnod[i][j] = k the j-th adjacent node to i-th node has number k
   JB
**/
void fixnodesel::create_subdom_graph_3d ()
{
  long i,j,k,l,m;
  long ned,nned,auxnedges;
  long *auxnned,*dupledges,*aux;
  long **auxedgenodes;
  double ts,te;
  
  ts  = clock();
  // number of all edges in mesh
  auxnedges = 0;
  // pomocne pole
  aux = new long[3];
  for(i = 0; i < ne; i++){
    // number of edge on i-th element - gtopology
    ned = top->give_ned(i);
    nned = top->give_nned(i);
    for(j = 0; j < ned; j++ ){
      top->give_edge_nodes(i,j,aux);
      m = 0;
      for(l = 0; l < nned; l++){
	if(nodmultip[aux[l]] > 1){
	  m++;
	}
      }
      if(m > 1){
	auxnedges++;
      }
    }
  }
  //fprintf(out,"pocet hran %ld\n",auxnedges);
  te  = clock();
  //fprintf (out,"pocet hran %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  // plneni pole s uzky na hranach
  // plneni s pole s poctem uzlu na hrane
  ts  = clock();
  auxnned = new long[auxnedges];
  auxedgenodes = new long*[auxnedges];
  dupledges = new long[auxnedges];
  k = 0;
  for(i = 0; i < ne; i++){
    ned = top->give_ned(i);
    nned = top->give_nned(i);
    for(j = 0; j < ned; j++ ){
      top->give_edge_nodes(i,j,aux);
      m = 0;
      for(l = 0; l < nned; l++){
	if(nodmultip[aux[l]] > 1){
	  m++;
	}
      }
      if(m > 1){
	auxnned[k] = nned;
	auxedgenodes[k] = new long[nned];
	for(l = 0; l < nned; l++){
	  auxedgenodes[k][l] = aux[l];
	}
	// trideni pole s uzly na hrane(od nejmensiho po nejvetsi)
	if(nned == 2){
	  if(auxedgenodes[k][0] > auxedgenodes[k][1]){
	    l = auxedgenodes[k][1];
	    auxedgenodes[k][1] = auxedgenodes[k][0];
	    auxedgenodes[k][0] = l;
	  }
	}
	if(nned == 3){
	  if(auxedgenodes[k][0] > auxedgenodes[k][2]){
	    l = auxedgenodes[k][2];
	    auxedgenodes[k][2] = auxedgenodes[k][0];
	    auxedgenodes[k][0] = l;
	  }
	  if(auxedgenodes[k][0] > auxedgenodes[k][1]){
	    l = auxedgenodes[k][1];
	    auxedgenodes[k][1] = auxedgenodes[k][0];
	    auxedgenodes[k][0] = l;
	  }
	  if(auxedgenodes[k][1] > auxedgenodes[k][2]){
	    l = auxedgenodes[k][2];
	    auxedgenodes[k][2] = auxedgenodes[k][1];
	    auxedgenodes[k][1] = l;
	  }
	}
	  dupledges[k] = 0;
	  k++;
      }
    }
  }
  delete []aux;
  //auxnedges = k;
  //fprintf(out,"pocet hran %ld\n",auxnedges);
  te  = clock();
  //fprintf (out,"seznam uzlu %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  
  
  // kontrolni tisk
  //for(i = 0; i < auxnedges; i++){
  //fprintf(out,"hrana %ld ma uzly:",i);
  //for(j = 0; j < auxnned[i]; j++ ){
  //fprintf(out,"   %ld",auxedgenodes[i][j]);
  //}
  //fprintf(out,"\n");
  //}
  ts  = clock();
  // hledani duplikovanych hran na hranici
  for(i = 0; i < auxnedges; i++){
    for(j = i+1; j < auxnedges; j++){
	if(auxedgenodes[i][0] == auxedgenodes[j][0]){
	  if(auxnned[i] == 2 && auxnned[j] == 2){
	    if(auxedgenodes[i][1] == auxedgenodes[j][1]){
	      dupledges[j]++;
	      //dupledges[i]++;
	    }
	  }
	  if(auxnned[i] == 3 && auxnned[i] == 3){
	    if(auxedgenodes[i][1] == auxedgenodes[j][1]){
	      if(auxedgenodes[i][2] == auxedgenodes[j][2]){
		dupledges[j]++;
		//dupledges[i]++;
	      }
	    }
	  }
	}
    }
  }
  te  = clock();
  //fprintf (out,"hledani dulikovanych hran %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  
  // kontrolni tisk
  //for(i = 0; i < auxnedges; i++){
  //fprintf(out,"%ld dupledges %ld\n",i,dupledges[i]);
  //}
  k = 0;
  for(i = 0; i < auxnedges; i++){
    if(dupledges[i] == 0){
      k++;
    }
  }
  //fprintf(out," k je %ld\n",k);
  ts  = clock();
    
  nadjacboundnod = new long[nbn];
  for(i = 0; i < nbn; i++){
    nadjacboundnod[i] = 0;
    for(j = 0; j < auxnedges; j++){
      if(dupledges[j] == 0) {
	for(k = 0; k < auxnned[j]; k++){
	  if(auxedgenodes[j][k] == lnbn[i]){
	    auxedgenodes[j][k] = i;
	    nadjacboundnod[i]++;
	  }
	}
      }
    } 
  }
  te  = clock();
  //fprintf (out,"pocet sousedu %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  
  //for(i = 0; i < nbn; i++){
  //fprintf(out,"%ld\n",nadjacboundnod[i]);
  //}  
  ts  = clock();
    
  adjacboundnod = new long*[nbn];
  for(i = 0; i < nbn; i++){
    adjacboundnod[i] = new long[nadjacboundnod[i]];
    m = 0;
    for(j = 0; j < auxnedges; j++){
      if(dupledges[j] == 0) {
	if(auxnned[j] == 2){
	  for(k = 0; k < auxnned[j]; k++){
	    if(auxedgenodes[j][k] == i){
	      if(k == 0){
		adjacboundnod[i][m] = auxedgenodes[j][1];
		m++;
	      }
	      if(k == 1){
		adjacboundnod[i][m] = auxedgenodes[j][0];
		m++;
	      }
	    }
	  }
	}
	if(auxnned[j] == 3){
	  for(k = 0; k < auxnned[j]; k++){
	    if(auxedgenodes[j][k] == lnbn[i]){
	      if(k == 0){
		adjacboundnod[i][m] = auxedgenodes[j][1];
		m++;
	      }
	      if(k == 1){
		adjacboundnod[i][m] = auxedgenodes[j][0];
		m++;
		adjacboundnod[i][m] = auxedgenodes[j][2];
		m++;
	      }
	      if(k == 2){
		adjacboundnod[i][m] = auxedgenodes[j][1];
		m++;
	      }
	    }
	  }
	}
      }
    }
  }
  // kontrolni tisk
  //for(i = 0; i < nbn; i++){
  // fprintf(out,"uzel %ld ma  %ld sousedu: ",i,nadjacboundnod[i]);
  //for(j = 0; j < nadjacboundnod[i]; j++){
  //fprintf(out,"   %ld",adjacboundnod[i][j]);
  //}
  //fprintf(out,"\n");
  //}
  //for(i = 0; i < nbn; i++){
  // fprintf(out,"uzel %ld ma lgnbn %ld\n",i,lgnbn[i]);
  //}
  te  = clock();
  //fprintf (out,"seznam sousedu %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  
  
  ts  = clock();
    
  //sorting 
  for(i = 0; i < nbn; i++){
    if(nadjacboundnod[i] == 2){
      if(adjacboundnod[i][0] > adjacboundnod[i][1]){
	j = adjacboundnod[i][1];
	adjacboundnod[i][1] = adjacboundnod[i][0];
	adjacboundnod[i][0] = j;
      }
    }
    else{
      // bubble sorting
      for(j = nadjacboundnod[i] - 1; j > 0; j--){
	for(k = 0; k < j; k++){
	  if(adjacboundnod[i][k] > adjacboundnod[i][k+1]){
	    m = adjacboundnod[i][k];   
	    adjacboundnod[i][k] = adjacboundnod[i][k+1];
	    adjacboundnod[i][k+1] = m;
	  }
	} 
      }
    }
  }
  te  = clock();
  //fprintf (out,"sorting %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  
  
  // kontrolni tisk
  // for(i = 0; i < nbn; i++){
  //     fprintf(out,"uzel %ld ma  %ld sousedu: ",i,nadjacboundnod[i]);
  //     for(j = 0; j < nadjacboundnod[i]; j++){
  //       fprintf(out,"   %ld",adjacboundnod[i][j]);
  //     }
  //     fprintf(out,"\n");
  //   }
  ts  = clock();
    
  // rewrite for creation of master graph
  for(i = 0; i < nbn; i++){
    for(j = 0; j < nadjacboundnod[i]; j++){
      k = adjacboundnod[i][j];
      adjacboundnod[i][j] = lgnbn[k];
    }
  }
  te  = clock();
  //fprintf (out,"prepis %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  
  fflush(out);
  
  for(i = 0; i < auxnedges; i++){
    delete []auxedgenodes[i];
  }
  delete []auxnned;
  delete []dupledges;
  delete []auxedgenodes;  
}

/**
   Function select minimal number of fixing nodes for 3D meshes
   Selection is based on graph theory and is described in dissertation JB
   @output coarseidentif - array with identification of boundaty nodes
   coarseidentif = 2 - remaining interface node
   coarseidentif = 3 - fixing node
   coarseidentif has tnbn components
   JB
 **/
void fixnodesel::select_minimal_number_3d()
{
  long i,j,a,b,c,d;
  long auxnv;
  
  if(myrank == 0){
    // identifikace vnitrnich uzlu na domene
    if (coarseidentif != NULL)
      delete []coarseidentif;
    
    coarseidentif = new long [tnbn];
    auxnv = 0;
    for(i = 0; i < tnbn; i++){
      coarseidentif[i] = 2;
      if(multip[i] > 2){
	//fprintf(out,"uzel %ld multip %ld\n",i,multip[i]);
	b = 0;
	c = 0;
	d = 0;
	for(j = 0; j < coarsenadjac[i]; j++){
	  a = coarseadjac[i][j];
	  //fprintf(out,"multiplicita %ld je %ld\n",coarseadjac[i][j],multip[a]);
	  if(multip[a] <  multip[i]){
	    b++;
	  }
	  if(multip[a] == multip[i]){
	    c++;
	  }
	  if(multip[a] > multip[i]){
	    d++;
	  }
	}
	//fprintf(out,"b %ld c %ld d %ld\n",b,c,d);
	// bod dotyku nodmultip[lnbn[i]] podoblasti - krizeni krivek o multiplicite > 2
	if(c == 0 && d == 0 && b == coarsenadjac[i]){
	  coarseidentif[i] = 3;
	}
	// pocatek krivky o multiplicite multip[i] 
	if(c == 1 && d == 0 ){
	  coarseidentif[i] = 3;
	}
      }
      if(coarseidentif[i] == 3){
	auxnv++;
      }
    }
    if(mespr == 1) fprintf(out,"Number of potencial fixing nodes in whole problem is %ld\n\n",auxnv);
    fflush(out);
  }
  
}

/**
   Function select fixing nodes for 3D meshes based on statistics only
   JB
 **/

void fixnodesel::select_fixing_nodes_on_master()
{
  long i,j,k;
  long *auxnv,*buff;
  long acticmultip,stop,buffsize,estMinCornod;
  MPI_Status stat;
  
  
  buffsize = maxnbn;
  buff = new long[buffsize];
  
  if(myrank == 0){
    estMinCornod = nproc*3;
    
    auxnv = new long[nproc];
    for(i = 0; i < nproc; i++){
      auxnv[i]=0;
    }
    // indentifikace na masteru 
    coarseidentif = new long[tnbn];
    // vyber uzlu s nejvetsi icmultiplicitou
    for(i = 0; i < tnbn; i++){
      if(multip[i] == maxmultip){
	coarseidentif[i] = 3;
	for(j = 0; j < multip[i]; j++){
	  auxnv[glinkdom[i][j]]++;
	}
      }
      else{
	coarseidentif[i] = 2;
      }
    }
    
    //for(i = 0; i < nproc; i++){
    //fprintf(out,"auxnv[%ld] = %ld\n",i,auxnv[i]);
    //}

    
    for(i = 0; i < nproc; i++){
      if(auxnv[i] < 3){
	stop = 0;
	acticmultip = maxmultip-1;
	while (stop != 1){
	  if(statdom[i][acticmultip-1] != 0){
	    for(j = 0; j < nbnd[i]; j++){
	      if(multip[cnbn[i][j]] == acticmultip){
		coarseidentif[cnbn[i][j]] = 3;
		for(k = 0; k < multip[cnbn[i][j]]; k++){
		  auxnv[glinkdom[cnbn[i][j]][k]]++;
		}
	      }
	    }
	    if(auxnv[i] >= 3){
		  stop = 1;
	    }
	    else{
	      stop = 0;
	      acticmultip--;
	    }
	  }
	  else{
	    acticmultip--;
	  }
	}
      }
    }
    
    // for(i = 0; i < nproc; i++){
    //       fprintf(out,"auxnv[%ld] = %ld\n",i,auxnv[i]);
    //     }
    
    j = 0;
    k = 0;
    for(i = 0; i < tnbn; i++){
      if(coarseidentif[i] == 3){
	j++;
      }
      else{
	k++;
      }
    }
    if(mespr == 1){
      fprintf(out,"\nNumber of fixing nodes in whole problem is %ld\n",j);
      fprintf(out,"\nNumber of boundary nodes in whole problem is %ld\n",k);
    }
    for(j = 1; j < nproc; j++){
      for(k = 0; k < nbnd[j]; k++){
	buff[k] = coarseidentif[cnbn[j][k]];
      }
      MPI_Send (buff,buffsize,MPI_LONG,j,myrank,MPI_COMM_WORLD);
    }
    
    nodeidentif = new long[nbn];
    for(j = 0; j < nbnd[0]; j++){
      nodeidentif[j] = coarseidentif[cnbn[0][j]];;
    }
    
    delete []auxnv;
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    
    nodeidentif = new long[nbn];
    for(i = 0; i < nbn; i++){
      nodeidentif[i] = buff[i];
    }
  }
  if(mespr == 1){
    for (i = 0; i < nbn; i++){
      if(nodeidentif[i] == 3) fprintf(out,"\n node number %ld is fixing node",lnbn[i]+1);
    }
    fflush(out);
    fprintf(out,"\n\n");
  }
  delete []buff;
}

/**
   Function selects fixing nodes based on geometrical condition of 3D mesh
   JB
 **/
void fixnodesel::select_fixing_nodes_geom_3d(long auxnv)
{
  long i,j,k,max;
  long *ver,first,second,third;
  double maxradius,radius,area,maxarea,cosalpha;
  double *contrad;
  long *vcomb;

  //double t1= clock();
  if(auxnv == 0){
    // inicialization of pseudorandom numbers
        srand((unsigned int) time(NULL));
        j = rand()%(nbn-1);
    
   //  max = 0;
//     for(i = 0; i < nbn; i++){
//       if(max < lgnbn[i]){
//     	max = lgnbn[i];
//     	j = i;
//       }
//     }
    //fprintf(out,"max je %ld\n\n",max);
    //fprintf(out,"prvni je %ld\n\n",lnbn[j]);
    
    vcomb = new long[3];
    vcomb[0] = lnbn[j];
    maxradius = 0.0;
    for(i = 0; i < nbn; i++){
      if(i != j){
	vcomb[1] = lnbn[i];
	radius = compute_length_of_vector(vcomb[0],vcomb[1]);
	//fprintf(out,"radius je %lf\n\n",radius);
	if(radius > maxradius){
	  maxradius = radius;
	  first = i;
	}
      }
    }
    nodeidentif[first]=3;
    //fprintf(out,"novy je %ld\n\n",lnbn[first]);
    //fprintf(out,"maxradius je %lf\n\n",maxradius);
    
    vcomb[0] = lnbn[first];
    maxradius = 0.0;
    for(i = 0; i < nbn; i++){
      if(i != first){
	vcomb[1] = lnbn[i];
	radius = compute_length_of_vector(vcomb[0],vcomb[1]);
	//fprintf(out,"radius je %lf\n\n",radius);
	if(radius > maxradius){
	  second = i;
	  maxradius = radius;
	}
      }
    }
    nodeidentif[second]=3;
    vcomb[0] = lnbn[first];
    vcomb[1] = lnbn[second];
    //if(mespr == 1) fprintf(out,"potencialni vrcholy jsou %ld %ld\n\n",lnbn[second],lnbn[first]);
    //fprintf(out,"maxradius je %lf\n\n",maxradius);
    
    
    maxarea = 0.0;
    for(i = 0; i < nbn; i++){
      if(vcomb[0] != lnbn[i] && vcomb[1] != lnbn[i]){
	vcomb[2] = lnbn[i];
	//fprintf(out,"treti je %ld\n\n",lnbn[i]);
	area = compute_area_of_triangle(vcomb[0],vcomb[1],vcomb[2]);
	//fprintf(out,"uzel %ld area je %lf\n\n",lnbn[i],area);
	if(maxarea < area){
	  maxarea = area;
	  third = i;
	}
      }
    }
    nodeidentif[third] = 3;
    //fprintf(out,"treti je %ld\n\n",lnbn[j]);
    //fprintf(out,"maxrarea je %lf\n\n",maxarea);
    delete []vcomb;
    
    //     for(i = 0; i < nbn; i++){
    //       if(nodeidentif[i] == 3){
    //     	fprintf(out,"%ld je fixing\n\n",lnbn[i]);
    //       }
    //     }
    //double t2=  clock();
    //fprintf(out,"cas hledani je %lf\n\n",(t2-t1)/(double)CLOCKS_PER_SEC);    
  }
  
  if(auxnv == 1){
    max = 0;
    for(i = 0; i < nbn; i++){
      if(nodeidentif[i] == 3){
	first = i;
      }
    }
    //fprintf(out,"max je %ld\n\n",max);
    //fprintf(out,"prvni je %ld\n\n",lnbndom[first]);
    
    vcomb = new long[3];
    print_err("array ver is not initialized nor even allocated!", __FILE__, __LINE__, __func__);    
    vcomb[0] = lnbn[ver[0]]; // ver is not even allocated nor initialized !!!???
    maxradius = 0.0;
    for(i = 0; i < nbn; i++){
      vcomb[1] = lnbn[i];
      radius = compute_length_of_vector(vcomb[0],vcomb[1]);
      //fprintf(out,"radius je %lf\n\n",radius);
      if(radius >= maxradius){
	  maxradius = radius;
	  second = i;
      }
    }
    //fprintf(out,"dalsi je %ld\n\n",lnbndom[ver[1]]);
    //fprintf(out,"maxradius je %lf\n\n",maxradius);
    //fprintf(out,"potencialni vrcholy jsou %ld %ld\n\n",lnbndom[second],lnbndom[first]);
    
    nodeidentif[second]=3;
    vcomb[0] = lnbn[first];
    vcomb[1] = lnbn[second];
    
    delete []ver;
    
    maxarea = 0.0;
    for(i = 0; i < nbn; i++){
      vcomb[2] = lnbn[i];
      area = compute_area_of_triangle(vcomb[0],vcomb[1],vcomb[2]);
      if(maxarea < area){
	maxarea = area;
      }
    }
    for(i = 0; i < nbn; i++){
      vcomb[2] = lnbn[i];
      area = compute_area_of_triangle(vcomb[0],vcomb[1],vcomb[2]);
      if(maxarea < area){
	j = i;
      }
    }
    nodeidentif[j] = 3;
    delete []vcomb;
    
    // 	for(i = 0; i < nbn; i++){
    // 	  if(nodeidentif[i] == 3){
    // 	    fprintf(out,"%ld je fixing\n\n",lnbn[i]);
    // 	  }
    // 	}
  }
  if(auxnv == 2){
    ver = new long[2];
    vcomb = new long[3];
    j = 0;
    k = 0;
    for(i = 0; i < nbn; i++){
      if(nodeidentif[i] == 3){
	ver[j] = lnbn[i];
	j++;
      }
      if(nodeidentif[i] == 2){
	k++;
      }
    }
    radius = compute_length_of_vector(ver[0],ver[1]);
    //fprintf(out,"radius %lf mezi uzly %ld %ld\n\n",radius,ver[0],ver[1]);
    maxradius = 0.0;
    j = 0;
    contrad = new double[k];
    vcomb[0]=ver[0];
    for(i = 0; i < nbn; i++){
      if(nodeidentif[i] == 2){
	vcomb[1] = lnbn[i];
	contrad[j] = compute_length_of_vector(vcomb[0],vcomb[1]);
	//fprintf(out,"uzel %ld contrad %lf\n\n",lnbndom[i],contrad[j]);
	if(maxradius < contrad[j]){
	  maxradius = contrad[j];
	}
	j++;
      }
    }
    //fprintf(out,"maxradius %lf\n\n",maxradius);
    j = 0;
    vcomb[0] = ver[1];
    for(i = 0; i < nbn; i++){
      if(nodeidentif[i] == 2){
	//fprintf(out,"uzel %ld contrad %lf maxradius %lf\n\n",lnbndom[i],contrad[j],maxradius);
	if(maxradius < contrad[j]){
	  vcomb[1] = lnbn[i];
	  radius = compute_length_of_vector(vcomb[0],vcomb[1]);
	  //fprintf(out,"radius %lf\n\n",radius);
	  if(radius > 0.05*sizefictdom){
	    vcomb[2] = ver[0];
	    cosalpha = compute_angle_of_vector_a_and_b(vcomb[0],vcomb[1],vcomb[2]);
	    if(cosalpha > 0.034906585 || cosalpha < 3.1066861){
	      nodeidentif[i] = 3;
	      //fprintf(out,"uzel  %ld bude oznacen jako fixing\n\n",lnbndom[i]);
	    }
	  }
	}
	j++;
	}
    }
    delete []contrad;
    delete []ver;
    delete []vcomb;
  }
  auxnv = 3;
}

/**
   Function check minimal number of selected fixing nodes for 3D meshes
   JB
 **/
void fixnodesel::check_minimal_number_3d()
{
  long a,b,i,j,k,l,m,auxnv;
  long buffsize,*buff;
  MPI_Status stat;
  long *nvdom;
  long **stor;
  
  buffsize = maxnbn+1;
  buff = new long[buffsize];
  if(myrank == 0){
    nvdom = new long[nproc];
    for(i = 0; i < nproc; i++){
      nvdom[i] = 0;
    }
    
    for(i = 0; i < tnbn; i++){
      if(coarseidentif[i] == 3){
	for(j = 0; j < multip[i]; j++){
	  k = glinkdom[i][j];
	  nvdom[k]++;
	}
      }
    }
    //     for(i = 0; i < nproc; i++){
    //       fprintf(out,"domena %ld nvdom %ld\n\n",i,nvdom[i]);    
    //     }
    
    //slave
    for(i = 1; i < nproc; i++){
      for(j = 0; j < nbnd[i]; j++){
	k = cnbn[i][j];
	buff[j] = coarseidentif[k];
      }
      buff[maxnbn] = nvdom[i];
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    //master
    for(j = 0; j < nbnd[0]; j++){
      k = cnbn[0][j];
      buff[j] = coarseidentif[k];
    }
    buff[maxnbn] = nvdom[0];
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  nodeidentif = new long[nbn];
  
  for(i = 0; i < nbn; i++){
    nodeidentif[i] = buff[i];
  }
  
  if(mespr == 1) fprintf(out,"Number of potential fixing nodes on subdomain is %ld\n\n",buff[maxnbn]);
  fflush(out);
  
  
  if(buff[maxnbn] < 3){
    auxnv = buff[maxnbn];
    select_fixing_nodes_geom_3d(auxnv);
    buff[maxnbn] = -3;
  }
  else{
    if(condfixing == nocondconer){
      check_triangle();
      auxnv = 0;
      for(i = 0; i < nbn; i++){
	if(nodeidentif[i] == 3){
	  auxnv++;
	}
      }
      if(auxnv != buff[maxnbn]){
	buff[maxnbn] = -1*auxnv;
      }
      if(auxnv < 3){
	select_fixing_nodes_geom_3d(auxnv);
	buff[maxnbn] = -3;
      }
    }
  }
  
  if(buff[maxnbn] < 0){
    for(i = 0; i < nbn; i++){
      buff[i] = nodeidentif[i]; 
    }
  }
  
  if(mespr == 1) fprintf(out,"Number of potential fixing nodes %ld on subdomain %d\n\n",abs(buff[maxnbn]),myrank);
  
  
  
  if(myrank == 0){
    //master
    stor = new long*[nproc];
    //fprintf(out,"domena 0 nvdom stare %ld nove %ld\n\n",nvdom[0],buff[maxnbn]);
    if (buff[maxnbn] < 0 ){
      nvdom[0] = -1*buff[maxnbn];
      stor[0] = new long [nbnd[0]];
      for(j = 0; j < nbnd[0]; j++){
      	stor[0][j] = buff[j];
      }
    }
    else{
      stor[0]= new long [nbnd[0]];
      for(j = 0; j < nbnd[0]; j++){
      	stor[0][j] = buff[j];
      }
    }
    
    //slave
    for(i = 1; i < nproc; i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);  
      k=domproc[stat.MPI_TAG];
      //fprintf(out,"domena %ld nvdom stare %ld nove %ld\n\n",k,nvdom[k],buff[maxnbn]);
      if (buff[maxnbn] < 0 ){
	nvdom[k] = -1*buff[maxnbn];
	stor[k] = new long [nbnd[k]];
	for(j = 0; j < nbnd[k]; j++){
	  stor[k][j] = buff[j];
	}
      }
      else{
	stor[k] = new long [nbnd[k]];
	for(j = 0; j < nbnd[k]; j++){
	  stor[k][j] = buff[j];
	}
      }
    }
    
    
    for(i = 0; i < tnbn; i++){
      // controlling
      m = 0;
      l = 0;
      if (coarseidentif[i] == 3){
	a = 1;
	b = 0;
      }
      else{
	a = 0;
	b = 1;
      }	
      for(j = 0; j < multip[i]; j++){
	if(coarseidentif[i] != stor[glinkdom[i][j]][glinknod[i][j]]){
	  if(stor[glinkdom[i][j]][glinknod[i][j]] == 3){
	    a++;
	  }
	  if(stor[glinkdom[i][j]][glinknod[i][j]] == 2){
	    b++;
	  }
	  m++;
	}
	if(nvdom[glinkdom[i][j]] <= 3){
	  l++;
	}
      }
      if( m != 0){
	//fprintf(out,"m = %ld a= %ld b = %ld\n",m,a,b);
	if(a >= b){
	  //fprintf(out,"prepis %ld z 2 na 3\n",i);
	  coarseidentif[i] = 3;
	  for(j = 0; j < multip[i]; j++){
	    if(stor[glinkdom[i][j]][glinknod[i][j]] == 2){
	      nvdom[glinkdom[i][j]]++;
	    }
	  }
	}
	else{
	  if(l == 0){
	    //fprintf(out,"prepis %ld z 3 na 2\n",i);
	    coarseidentif[i] = 2;
	    for(j = 0; j < multip[i]; j++){
	      if(stor[glinkdom[i][j]][glinknod[i][j]] == 3){
		nvdom[glinkdom[i][j]]--;
	      }
	    }
	  }
	  else{
	    
	  }
	}
      }
    }
    if(mespr == 1){
      for(i = 0; i < nproc; i++){
	fprintf(out,"\n\n\n Number of fixing nodes on domain %ld is %ld\n",i,nvdom[i]);
      }
    }
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);  

  delete []buff;

  if(mespr == 1){
    for(i = 0; i < nbn; i++){
      if(nodeidentif[i] == 3){
	fprintf(out,"local number of fixng node: %ld\n\n",lnbn[i]);
      }
    }
  }
    


  //info
  if(myrank == 0){
    m = 0;
    for(i = 0; i < nproc; i++){
      nvdom[i] = 0;
    }
    for(i = 0; i < tnbn; i++){
      if(coarseidentif[i] == 3){
	m++;
	for(j = 0; j < multip[i]; j++){
	  k = glinkdom[i][j];
	  nvdom[k]++;
	}
      }
    }
    if(mespr == 1){
      fprintf(out,"\n\n\n Total number of fixing nodes in whole problem is %ld\n",m);
      for(i = 0; i < nproc; i++){
	fprintf(out,"\n\n\n Number of fixing nodes on domain %ld is %ld\n",i,nvdom[i]);
      }
    }
    
    delete []nvdom;
    for(i = 0; i < nproc; i++){
      if(stor[i] != NULL){
	delete []stor[i];
      }
    }
    delete []stor;
  }
}

/**
   Function creates boundary graph on master processor
   @output coarsenadjac - number of boundary nodes adjacent to boundary node
   coarsenadjac has tnbn components
   coarsenadjac[i] = j - the i-th global node has j adjacent nodes
   @output coarseadjac - list of boundaru nodes adjacent to boundary node
   coarseadjac has tnbn components
   coarseadjac[i] has multip[i] components
   coarseadjac[i][j] = k - the j-th adjacent node to i-th node has global number k
   JB
 **/
void fixnodesel::create_master_graph()
{
  long i,j,k,m,a,b;
  long buffsize,max,adr1,adr2;
  long *buff,**aux,*naux,**auxnadjac,*help;
  MPI_Status stat;
  
  buffsize = nbn;
  for(i = 0; i < nbn; i++){
    buffsize += nadjacboundnod[i];
  }
  
  //fprintf(out,"buffsize = %ld\n",buffsize);

  
  if (myrank==0){
    //  master contribution
    max = buffsize;
    
    for (i=1;i<nproc;i++){
      MPI_Recv (&k,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //  slave contributions
      if (max<k){
	max=k;
      }
    }
    
    buffsize = max;
    
    for (i=1;i<nproc;i++){
      MPI_Send (&buffsize,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  // ********
  // slaves
  // ********
  else{
    MPI_Send (&buffsize,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (&buffsize,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  //fprintf(out,"buffsize = %ld\n",buffsize);
  
  buff = new long[buffsize];
  
  buff[0] = 0;
  for(i = 1; i < nbn; i++){
    buff[i] = buff[i-1]+nadjacboundnod[i-1];
  }
  
  m = nbn;
  for(i = 0; i < nbn; i++){
    for(j = 0; j < nadjacboundnod[i]; j++){
      buff[m]  = adjacboundnod[i][j];
      m++;
    }
  }
  for(i = m; i < buffsize; i++){
    buff[m] = -1;
    m++;
  }
  
  //for(i = 0; i < buffsize+1; i++){
  //fprintf(out,"kont buff %ld %ld\n",i,buff[i]);
  //}
  
  //delete nadjacboundnod
  // opravdu je to nutne?, nebudu to potrebovat?
  for(i = 0; i < nbn; i++){
    delete []adjacboundnod[i];
  }
  delete []adjacboundnod;
  delete []nadjacboundnod;
  
  if (myrank == 0){
    auxnadjac = new long*[nproc];
    aux = new long*[nproc];
    naux = new long[nproc];
    // master contribution
    k=domproc[0];
    auxnadjac[k] = new long[nbnd[k]+1];
    for(j = 0; j < nbnd[k]; j++){
      auxnadjac[k][j] = buff[j];
    }
    naux[k] = 0;
    m = 0;
    for(j = nbnd[k]; j < buffsize; j++){
      if(buff[j] != -1){
	naux[k]++;
      }
      else{
	m++;
      }
      if(m == 2){
	break;
      }
    }
    //fprintf(out,"%ld naux = %ld %ld\n",k,naux[k],j);
    aux[k] = new long[naux[k]];
    auxnadjac[k][nbnd[k]] = naux[k];
    for(j = 0; j < naux[k]; j++){
      aux[k][j] = buff[j+nbnd[k]];
    }
    
    
    //slaves contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=domproc[stat.MPI_TAG];
      auxnadjac[k] = new long[nbnd[k]+1];
      for(j = 0; j < nbnd[k]; j++){
	auxnadjac[k][j] = buff[j];
      }
      naux[k] = 0;
      m = 0;
      for(j = nbnd[k]; j < buffsize; j++){
	if(buff[j] != -1){
	  naux[k]++;
	}
	else{
	  m++;
	}
	if(m == 2){
	  break;
	}
      }
      //fprintf(out,"%ld naux = %ld %ld\n",k,naux[k],j);
      aux[k] = new long[naux[k]];
      auxnadjac[k][nbnd[k]] = naux[k];
      for(j = 0; j < naux[k]; j++){
	aux[k][j] = buff[j+nbnd[k]];
      }
    }
    
    //for(i = 0; i < nproc; i++){
    //fprintf(out,"domain %ld %ld:",i,naux[i]);
    //for(j = 0; j < nbnd[i]; j++){
    //fprintf(out,"   %ld",auxnadjac[i][j]);
    //}
    //fprintf(out,"\n");
    //}
    
    //for(i = 0; i < nproc; i++){
    //fprintf(out,"domain %ld\n",i);
    //for(j = 0; j < nbnd[i]; j++){
    //adr1=auxnadjac[i][j];
    //adr2=auxnadjac[i][j+1];
    //fprintf(out,"uzel %ld pocet sousedu: %ld ",j,adr2-adr1);
    //for(k = adr1; k < adr2; k++){
    //  fprintf(out," %ld",aux[i][k]);
    //}
    // fprintf(out,"\n");
    //}
    //}
    
    max = 0;
    for(i = 0; i < tnbn; i++){
      b = 0;
      for(j = 0; j < multip[i]; j++){
      adr1=auxnadjac[glinkdom[i][j]][glinknod[i][j]];
      adr2=auxnadjac[glinkdom[i][j]][glinknod[i][j]+1];
      a = adr2-adr1;
      b+=a;
      }
      if(max < b){
	max = b;
      }
    }
    //fprintf(out,"max je %ld\n",max);    
    
    coarsenadjac = new long[tnbn];
    coarseadjac = new long*[tnbn];
    help = new long[max];
    for(i = 0; i < tnbn; i++){
      //fprintf(out,"uzel  %ld\n",i);    
      b = 0;
      for(j = 0; j < multip[i]; j++){
	adr1=auxnadjac[glinkdom[i][j]][glinknod[i][j]];
	adr2=auxnadjac[glinkdom[i][j]][glinknod[i][j]+1];
	for(k = adr1; k < adr2; k++){
	  help[b] = aux[glinkdom[i][j]][k];
	  b++;
	}
      }
      //fprintf(out,"b je %ld\n\n",b);    
      //print
      //for(j = 0; j < b; j++){
      //fprintf(out,"    %ld",help[j]);    
      //}
      //fprintf(out,"\n"); 
      if(b == 2){
	if(help[0] != help[1]){
	  par_print_err(myrank,"problem with node %ld\n\n", __FILE__, __LINE__, __func__);
	}
	else{
	  coarsenadjac[i] = 1;
	  coarseadjac[i] = new long[coarsenadjac[i]];
	  coarseadjac[i][0] = help[0];
	  //fprintf(out,"ps %ld  s %ld\n",coarsenadjac[i],coarseadjac[i][0]); 
	}
      }
      else{
	// bubble sorting
	for(j = b - 1; j > 0; j--){
	  for(k = 0; k < j; k++){
	    if(help[k] > help[k+1]){
	      m = help[k];   
	      help[k] = help[k+1];
	      help[k+1] = m;
	    }
	  } 
	}
	// print
	//for(j = 0; j < b; j++){
	//fprintf(out,"    %ld",help[j]);    
	//}
	//fprintf(out,"\n"); 
	
	coarsenadjac[i] = 1;
	for(j = 1; j < b; j++){
	  if(help[j] != help[j-1]){
	    coarsenadjac[i]++;
	  }
	}
	//fprintf(out,"s %ld\n",coarsenadjac[i]); 
	coarseadjac[i] = new long[coarsenadjac[i]];
	coarseadjac[i][0] = help[0];
	m = 1;
	for(j = 1; j < b; j++){
	  if(help[j] != help[j-1]){
	    coarseadjac[i][m] = help[j];
	    m++;
	  }
	}
      }
    }
    delete []help;
    
//     for(i = 0; i < tnbn; i++){
//       fprintf(out,"uzel %ld ma  %ld sousedu:",i,coarsenadjac[i]); 
//       for(j = 0; j < coarsenadjac[i]; j++){
//         fprintf(out,"    %ld",coarseadjac[i][j]);    
//       }
//       fprintf(out,"\n");    
//     }
    
    for(i = 0; i < nproc; i++){
      delete []aux[i];
      delete []auxnadjac[i];
    }
    delete []aux;
    delete []naux;
    delete []auxnadjac;
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
}

/**
   Function creates new local to global correspondence related to the selected fixing nodes
   Function changes mesh desctiption to boundary_node
   @param ltg - array with local to global corresponcence
   @output rewritten ltg
   ltg = -1 - internal node
   ltg < -1 - remaining interface node
   ltg > 1 - fixing node
   JB
 **/
void fixnodesel::assemble_ltg_with_fixings(long *ltg)
{
  long i,k;
  long buffsize,c,b;
  long *buff,*coarsenumb;
  MPI_Status stat;
  
  
  
  buffsize = maxnbn;
  buff = new long[buffsize];
  
  //  numbering of fixing and edge nodes on master
  if(myrank==0){
     coarsenumb = new long[tnbn];
    // numbering of fixing and edge nodes on master
    c = 1;
    b = 1;
    for(i = 0; i < tnbn; i++){
      //fprintf(out,"%ld %ld\n",i+1,coarseidentif[i]);
      if(coarseidentif[i] == 3){
        coarsenumb[i] = c*(-1);
        c++;
	//fprintf(out,"now %ld %ld\n",i+1,coarsenumb[i]);
      }
      if(coarseidentif[i] == 2){
	coarsenumb[i] = b;
	b++;
	//fprintf(out,"now %ld %ld\n",i+1,coarsenumb[i]);
      }
    }
    if(mespr == 1){
      fprintf(out,"\n\n\n Number of fixing nodes in whole problem is %ld\n",c-1);
      fprintf(out,"\n\n\n Number of boundary nodes in whole problem is %ld\n",b-1);
    }
    if((c-1)+(b-1) != tnbn){
      par_print_err(myrank,"Problem: the number of fixing and remainning nodes is not equal to total number of boundary nodes\n", __FILE__, __LINE__, __func__);
      if(mespr == 1) fprintf (out,"Problem: the number of fixing and remainning nodes is not equal to total number of boundary nodes (file %s, line %d).\n",__FILE__,__LINE__);
      
    }
    // For Slaves
    for (i = 1;i < nproc; i++){
      for(k = 0; k < nbnd[i]; k++){
        buff[k] =  coarsenumb[cnbn[i][k]];
      }
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
    // For master
    for(k = 0; k < nbnd[0]; k++){
      buff[k] =  coarsenumb[cnbn[0][k]];
    }
    delete []coarsenumb;
  }
  
  // slave
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
//   fprintf(out,"\n\nfixing nodes - buff\n");
//    for(i = 0; i < nbn; i++){
//      fprintf(out,"%ld %ld\n",i+1,buff[i]);
//    }
   
  for(i = 0; i < nbn; i++){
    ltg[lnbn[i]]=buff[i]-1;
  }
  nin = nn - nbn;
  for(i = 0; i < nin; i++){
    ltg[lnin[i]]=-1;
  }
  
  // for(i = 0; i < nn; i++){
//     fprintf(out,"ltg[%ld] = %ld\n",i+1,ltg[i]);
//   }
  
  fflush(out);
  // deleting 
  delete []buff;
}


/**
   Function computes centers of gravity of boundary surfaces
   @output nsurfcentres - the number of centers of gravity of boundary surface
   nsurfcenters has nsurf components
   nsurcenters[i] = j - the i-th surface has j-th centres
   @output surfcenters - the list of centers of gravity of boundary surface
   surfcenters has nsurf componets
   surfcenters[i] has nsurfcenters components
   surfcenters[i][j]=k - the j-th centre of i-th surfaces has global number k
   JB
**/
void fixnodesel::select_centers_of_surfaces()
{
  long i,j,k,l,max;
  long buffsize,*buff,*aux;
  MPI_Status stat;
  double t1,t2,tc;

  long n;
  long **centres;
  long *ncentres;
  double x,y,z,cx,cy,cz,nc;
  double min;
  double *lengths;
  long **auxnodes;
  long subdomnsurf;
  long *subdomnsurfmem;
  long **subdomsurfmem;
  
  // maximal number of centres
  long maxncentres;
  maxncentres = 5;
  
  buffsize = maxnbn;
  buff = new long[buffsize];
  
  // surfnod if the nodes lie more than 2 subdomains than surfnod = -1, else surfnod = surafce number
  if(myrank == 0){
    for(i = 1; i < nproc; i++){
      for(j = 0; j < nbnd[i]; j++){
	buff[j] = surfnod[cnbn[i][j]];
      }
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
      
    }
    for(j = 0; j < nbnd[0]; j++){
      buff[j] = surfnod[cnbn[0][j]];
    }
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);      
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  max = 0;
  for(i = 0; i < nbn; i++){
    if(max < buff[i]){
      max = buff[i];
    }
  }
  //fprintf(out,"max je %ld\n",max);

  max++;
  aux = new long[max];
  for(i = 0; i < max; i++){
    aux[i] = 0;
  }
  
  for(i = 0; i < nbn; i++){
    if(buff[i] != -1){
      j = buff[i];
      aux[j]++;
    }
  }
  
  
  subdomnsurf = 0;
  for(i = 0; i < max; i++){
    if(aux[i] != 0){
      subdomnsurf++;
    }
  }
  
  subdomnsurfmem = new long[subdomnsurf];
  subdomsurfmem  = new long*[subdomnsurf];
  for(i = 0; i < subdomnsurf; i++){
    subdomnsurfmem[i] = 0;
  }
  
  j = 0;
  for(i = 0; i < max; i++){
    if(aux[i] != 0){
      subdomsurfmem[j] = new long[aux[i]];
      aux[i] = j;
      j++;
    }
  }
  
  for(i = 0; i < nbn; i++){
    if(buff[i] != -1){
      j = buff[i];
      k = aux[j];
      subdomsurfmem[k][subdomnsurfmem[k]]=i;
      subdomnsurfmem[k]++;
      buff[i] = -2;
    }
  }
  delete[]aux;
    
  
  // fprintf(out,"podoblast %ld ma %ld ploch\n",myrank,subdomnsurf);
  //   for(i = 0; i < subdomnsurf; i++){
  //     fprintf(out,"plocha %ld ma %ld uzlu\n",i,subdomnsurfmem[i]);
  //     for(j = 0; j < subdomnsurfmem[i]; j++){
  //       fprintf(out,"   %ld",lgnbn[subdomsurfmem[i][j]]);
  //     }
  //     fprintf(out,"\n");
  //   }
  
  max = 0;
  for(i = 0; i < subdomnsurf; i++){
    if(max < subdomnsurfmem[i]){
      max = subdomnsurfmem[i];
    }
  }
  
  //fprintf(out,"max je %ld\n",max);
  lengths = new double[max];
  ncentres = new long[subdomnsurf];
  centres = new long*[subdomnsurf];
  realcg = new double*[subdomnsurf];
  for(i = 0; i < subdomnsurf; i++){
    realcg[i] = new double[3];
    tc = 0.0;
    t1= clock();
    x = 0.0;
    y = 0.0;
    z = 0.0;
    nc = 0.0;
    for(j = 0; j < subdomnsurfmem[i]; j++){
      x += top->gnodes[lnbn[subdomsurfmem[i][j]]].x;
      y += top->gnodes[lnbn[subdomsurfmem[i][j]]].y;
      z += top->gnodes[lnbn[subdomsurfmem[i][j]]].z;
      nc += 1.00;
    }
    // centre of gravity
    //     cx=x/(double)subdomnsurfmem[i];
    //     cy=y/(double)subdomnsurfmem[i];
    //     cz=z/(double)subdomnsurfmem[i];
    cx=x/nc;
    cy=y/nc;
    cz=z/nc;
    realcg[i][0]=cx;
    realcg[i][1]=cy;
    realcg[i][2]=cz;
    // fprintf(out,"Coorodinate of the centre of %ld surface\n",i);
    //     fprintf(out,"cx %lf\n",cx);
    //     fprintf(out,"cy %lf\n",cy);
    //     fprintf(out,"cz %lf\n",cz);
    min=DBL_MAX;
    //x = top->gnodes[lnbn[subdomsurfmem[i][0]]].x;
    //y = top->gnodes[lnbn[subdomsurfmem[i][0]]].y;
    //z = top->gnodes[lnbn[subdomsurfmem[i][0]]].z;
    //min=(cx-x)*(cx-x)+(cy-y)*(cy-y)+(cz-z)*(cz-z);
    //fprintf(out,"min %lf\n",min);
    for(j = 0; j < subdomnsurfmem[i]; j++){
      x = top->gnodes[lnbn[subdomsurfmem[i][j]]].x;
      y = top->gnodes[lnbn[subdomsurfmem[i][j]]].y;
      z = top->gnodes[lnbn[subdomsurfmem[i][j]]].z;
      lengths[j]=(cx-x)*(cx-x)+(cy-y)*(cy-y)+(cz-z)*(cz-z);
      //lengths[j]=sqrt(lengths[j]);
      //fprintf(out,"lengths %.20le\n",lengths[j]);
      if(lengths[j] < min){
	min = lengths[j];
	//fprintf(out,"min %.20le\n",min);
      }
    }
    //fprintf(out,"minimal distance from centre to node %.20le\n",min);
    
    
    double diff;
    ncentres[i] = 0;
    for(j = 0; j < subdomnsurfmem[i]; j++){
      diff=fabs(min-lengths[j]);
      if(diff < tol){
 	ncentres[i]++;
      }
    }
    
    //fprintf(out,"ncenters[i] %ld\n",ncentres[i]);
    centres[i] = new long[ncentres[i]];
    n = 0;
    // centres are in global numbers
     for(j = 0; j < subdomnsurfmem[i]; j++){
       diff=fabs(min-lengths[j]);
       if(diff < tol){
	 //fprintf(out,"%ld %lf %lf %lf\n",lnbn[subdomsurfmem[i][j]],top->gnodes[lnbn[subdomsurfmem[i][j]]].x,top->gnodes[lnbn[subdomsurfmem[i][j]]].y,top->gnodes[lnbn[subdomsurfmem[i][j]]].z);
	 centres[i][n] = subdomsurfmem[i][j];//lgnbn[subdomsurfmem[i][j]];
	 n++;
       }
     }
     t2= clock();
     //   fprintf(out,"cas vypoctu centra je %lf\n\n",(t2-t1)/(double)CLOCKS_PER_SEC);    
     tc += ((t2-t1)/(double)CLOCKS_PER_SEC);
     //only for testing
     //ncentres[i] = 1;
     //centres[i] = new long[ncentres[i]];
     //centres[i][0] = subdomsurfmem[i][0];
  }
  delete []lengths;
 
  
  //fprintf(out,"cas vypoctu center je %lf\n\n",(tc)/(double)CLOCKS_PER_SEC);    
  // for(i = 0; i < nbn; i++){
  //     fprintf(out,"lgnbn   %ld\n",lgnbn[i]);    
  //   }
  
  if(mespr == 1){
    for(i = 0; i < subdomnsurf; i++){
      fprintf(out,"Surface %ld has %ld centres:\n",i,ncentres[i]);    
      for(j = 0; j < ncentres[i]; j++){
	fprintf(out,"   %ld",centres[i][j]);    
      }
      fprintf(out,"\n");    
    }
  }
  
  
  for(i = 0; i < subdomnsurf; i++){
    for(j = 0; j < ncentres[i]; j++){
      k = centres[i][j];
      buff[k] = 1;
    }
  }
  
//   for(i = 0; i < nbn; i++){
//     fprintf(out,"buff[%ld] %ld\n",i,buff[i]);    
//   }
  
  if(myrank == 0){
    auxnodes = new long*[nproc];
    auxnodes[0] = new long[nbnd[0]];
    for(j = 0; j < nbnd[0]; j++){
      auxnodes[0][j] = buff[j];
    }
    for(i = 1; i < nproc; i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k = domproc[stat.MPI_TAG];
      auxnodes[k] = new long[nbnd[k]];
      for(j = 0; j < nbnd[k]; j++){
	auxnodes[k][j] = buff[j];
      }
    }
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete []buff;
  
  if(myrank == 0){
    
    // for(i = 0; i < nproc; i++){
//       for(j = 0; j < nbnd[i]; j++){
// 	fprintf(out,"auxnodes %ld\n",auxnodes[i][j]);
//       }
//     }
    
    nsurfcentres = new long[nsurf];
    surfcenters = new long*[nsurf];
    for(i = 0; i < nsurf; i++){
      nsurfcentres[i] = 0;
    }
    
    for(i = 0; i < nproc; i++){
      for(j = 0; j < nbnd[i]; j++){
	if(auxnodes[i][j] == 1){
	  k = cnbn[i][j];
	  l = surfnod[k];
	  if(glinkdom[k][1] <= i){
	    nsurfcentres[l]++;
	  }
	}
      }
    }
    //fprintf(out,"nsurf %ld\n",nsurf);    
    
    for(i = 0; i < nsurf; i++){
      //fprintf(out,"Plocha ma %ld tezist %ld\n",i,nsurfcentres[i]);
      surfcenters[i] = new long[nsurfcentres[i]];
      nsurfcentres[i] = 0;
    }
    
    for(i = 0; i < nproc; i++){
      for(j = 0; j < nbnd[i]; j++){
	if(auxnodes[i][j] == 1){
	  k = cnbn[i][j];
	  l = surfnod[k];
	  if(glinkdom[k][1] <= i){
	    surfcenters[l][nsurfcentres[l]]=k;
	    nsurfcentres[l]++;
	  }
	}
      }
    }
    if(mespr == 1){
      for(i = 0; i < nsurf; i++){
	fprintf(out,"Surface %ld has %ld centres:\n",i,nsurfcentres[i]);
	for(j = 0; j < nsurfcentres[i]; j++){
	  fprintf(out,"   %ld",surfcenters[i][j]);
	}
	fprintf(out,"\n");
      }
    }
    for(i = 0; i < nproc; i++){
      delete []auxnodes[i];
    }
    delete []auxnodes;
  }
  
  for(i = 0; i < subdomnsurf; i++){
    delete []centres[i];
  }
  delete []centres;
  delete []ncentres;
  
}


/**
   Function marks nodes on boundary surfaces
   @output surfnodmark - list of marks of boundary nodes on boundary surfaces
   surfnodmark has nsurf components
   surfnodmark[i] has nsurfmember components
   surfnodmark[i][j] = k - the j-th member of i-th surface has mark k
   k is the path length form centre of gravity of boundary surface to j-th member
   JB
 **/
void fixnodesel::mark_surfnodes()
{
  long i,j,k,l,m,n,a,b;
  
  surfnodmark = new long*[nsurf];
  for(i = 0; i < nsurf; i++){
    //fprintf(out,"plocha %ld\n",i);
    surfnodmark[i] = new long[nsurfmembers[i]];
    
    for(j = 0; j < nsurfmembers[i]; j++){
      surfnodmark[i][j] = -1;
    }
    
    for(j = 0; j < nsurfcentres[i]; j++){
      k = surfcenters[i][j];
      for(l = 0; l < nsurfmembers[i]; l++){
	if(surfmembers[i][l] == k){
	  surfnodmark[i][l] = 0;
	  break;
	}
      }
    }
    
    m = 0;
    n = 0;
    for(m = 0; m < nsurfmembers[i]; m++){
      for(j = 0; j < nsurfmembers[i]; j++){
	if(surfnodmark[i][j] == m){
	  a = surfmembers[i][j];
	  b = coarsenadjac[a];
	  //fprintf(out,"uzel %ld sousedu: %ld\n",a,b);
	  for(k = 0; k < b; k++){
	    l = coarseadjac[a][k];
	    //fprintf(out,"soused %ld\n",l);
	    if(surfnod[l] == i && surfnodmark[i][surfnodpoint[l]] == -1){
	      surfnodmark[i][surfnodpoint[l]] = m+1;
	      n++;
	    }
	  }
	}
      }
    }
    
    // for(j = 0; j < nsurfmembers[i]; j++){
    //       fprintf(out,"%ld %ld\n",surfmembers[i][j],surfnodmark[i][j]);
    //     }
    
  }
}

void fixnodesel::mark_ring_nodes()
{
  
}

/**
   Function selects optimal number of fixing nodes for 2D meshes
   @output methodcondcor - place of fixing node condensation (curvecond)
   @output typecondcur - method of condensation of fixing nodes on boundary curves
   typecondcur = nth_member - nedges[i] is prime number
   typecondcur = n_part_curve - nedges[i] is non-prime number
   @output automember - list of nmembers for additon of fixing nodes
   automember has ncurves components
   automember[i] = j - i-th curve will be cut onto j-th parts in case of  n_part_curve
   or
   automember[i] = j - on i-th curve will be selected each j-th node in case of  nth_member
   
   The optimum number of fixing nodes was established as 5% of tnbn
   JB
 **/
void fixnodesel::select_optimum_2d()
{
  double optimum;
  double curopt;
  long nodes;
  long i, edges,tnedges;
  long stop,e,lnmember,ledges;//,lnedges;
  double zbytek,rat;
  bool prime;
  bool nondivide;
  
  
  tnedges = 0;
  for(i = 0; i < ncurves; i++){
    tnedges+=nedges[i];
  }
  //fprintf(out,"tnedges %ld\n",tnedges);
  methodcondcor = curvecond;
  typecondcur=n_part_curve;
  automember = new long[ncurves];
  // 0.05 - heuristic aproximation
  optimum = 0.05*tnbn;
  optimum = round_it(optimum);
  if(mespr == 1) fprintf(out,"Optimum is %ld nodes\n",(long)optimum);
  for(i = 0; i < ncurves; i++){
    nondivide = false;
    prime = isPrime(nedges[i]);
    switch(nondivide){
    case true:
      edges=nedges[i]+1;
      break;
    case false:
      edges=nedges[i];
      break;
    }
    rat=(double)nedges[i]/tnedges;
    curopt=rat*optimum;
    curopt = round_it(curopt);
    nodes = (long)curopt;
    //fprintf(out,"optimum for  curve %ld is %ld nodes\n",i,nodes);
    nodes++;
    if(edges%nodes == 0) {
      automember[i] = nodes;
    }
    else{
      lnmember = nodes;
      ledges = edges;
      // nova volba nmember
      stop = 0;
      e = 0;
      while (stop == 0){
	zbytek = ledges % lnmember;
	if( zbytek != 0){
	  //fprintf(out,"zbytek = %lf\n",zbytek);
	  if(e == 0){
	    lnmember--;
	  }
	  if(e == 1){
	    lnmember++;
	  }
	  //fprintf(out,"lnmember = %ld\n",lnmember);
	  stop = 0;
	}
	else{
	  //fprintf(out,"nova volba nmember = %ld\n",lnmember);
	  stop = 1;
	}
	if( ledges == lnmember){
	  e = 1;
	  lnmember = nmembercur;
	  stop = 0;
	}
	if(lnmember == 1){
	  stop = 1;
	  e = 2;
	  //fprintf(out,"%ld je prvocislo, nelze delit - nelze vybrat kazdy n-ty\nbude upraven pocet hran na edges+1 %ld\n",ledges);
	  ledges++;
	  lnmember = nmembercur;
	}
      }
      if(e != 2){
	automember[i] = lnmember;
      }
      else{
	automember[i] = 1;
      }
    }
    if(mespr == 1) fprintf(out,"curve %ld will be split into %ld parts\n",i,automember[i]);
  }

  
//   switch(nondivide){

//   case true:{
//     methodcondcor = curvecond;
//     typecondcur = nth_memb ;
//     automember = new long[ncurves];
//     // 0.05 - heuristic aproximation
//     optimum = 0.05*tnbn;
//     optimum = round_it(optimum);
//     fprintf(out,"Optimum is %ld nodes\n",(long)optimum);
//     for(i = 0; i < ncurves; i++){
//       curopt = (((double)nmembers[i]/(double)tnbn)*optimum);
//       curopt = round_it(curopt);
//       nodes = (long)curopt;
//       nodes =nodes+1;
//       automember[i]=(long)nmembers[i]/nodes;
//     }
//     for(i = 0; i < ncurves; i++){
//       fprintf(out,"curve %ld will be added node on each %ld member\n",i,automember[i]);
//     }
//     break;
//   }
//   case false:{
//     methodcondcor = curvecond;
//     typecondcur=n_part_curve;
//     automember = new long[ncurves];
//     // 0.05 - heuristic aproximation
//   optimum = 0.05*tnbn;
//   optimum = round_it(optimum);
//   fprintf(out,"Optimum is %ld nodes\n",(long)optimum);
//   for(i = 0; i < ncurves; i++){
//     curopt = (((double)nmembers[i]/(double)tnbn)*optimum);
//     curopt = round_it(curopt);
//     nodes = (long)curopt;
//     nodes =nodes+1;
//     fprintf(out,"optimum for  curve %ld is %ld nodes\n",i,nodes);
//     if(nedges[i]%nodes == 0) {
//       automember[i] = nodes;
//     }
//     else{
//       lnmember = nodes;
//       lnedges = nedges[i];
//       // nova volba nmember
//       stop = 0;
//       e = 0;
//       while (stop == 0){
// 	zbytek = lnedges % lnmember;
// 	if( zbytek != 0){
// 	  //fprintf(out,"zbytek = %lf\n",zbytek);
// 	  if(e == 0){
// 	    lnmember--;
// 	  }
// 	  if(e == 1){
// 	    lnmember++;
// 	  }
// 	  //fprintf(out,"lnmember = %ld\n",lnmember);
// 	  stop = 0;
// 	}
// 	else{
// 	  //fprintf(out,"nova volba nmember = %ld\n",lnmember);
// 	  stop = 1;
// 	}
// 	if( lnedges == lnmember){
// 	  e = 1;
// 	  lnmember = nmembercur;
// 	  stop = 0;
// 	}
// 	if(lnmember == 1){
// 	  stop = 1;
// 	  e = 2;
// 	  //fprintf(out,"%ld je prvocislo, nelze delit - nelze vybrat kazdy n-ty\nbude upraven pocet hran na nedges+1 %ld\n",lnedges);
// 	  lnedges++;
// 	  lnmember = nmembercur;
// 	}
//       }
//       if(e != 2){
// 	automember[i] = lnmember;
//       }
//       else{
// 	automember[i] = 1;
//       }
//     }
//   }
  
  
//   for(i = 0; i < ncurves; i++){
//     fprintf(out,"curve %ld will be split into %ld parts\n",i,automember[i]);
//   }
//   break;
//   }
//   }
  
}

/**
   The most important function of this @class
   Function manages selection of fixing nodes
   JB
 **/
void fixnodesel::fixing_detection(long *ltg)
{
  double ts,te,tc;
  tc = 0.0;
  ts  = clock();
  
  // set of fictitious subdomain
  set_fictitious_subdomain();
  te = clock();
  tc = tc+te-ts;
  if(mespr == 1) fprintf (out,"Time of set of fictitious subdomain %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  // size of fictitious subdomain
  ts  = clock();
  compute_size_of_fictitious_subdomain();
  te = clock();
  tc = tc+te-ts;
  if(mespr == 1) fprintf (out,"Time of computing of size of fictitious subdomain %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);
  // set of spatial dimension of problem
  ts  = clock();
  give_whole_dim();
  te = clock();
  tc = tc+te-ts;
  if(mespr == 1)  fprintf (out,"Time of set of spatial dimension %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
  fprintf(out,"\n\n\nDimension of problem is %ld\n",dim);
  ts = clock();
  create_link_dom_nod();
  te = clock();
  tc = tc+te-ts;
  if(mespr == 1)  fprintf (out,"Time of creation of globloc link %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
  switch(dim){
  case -1:
    par_print_err(myrank,"Problem:Spatial dimension is not detected\n", __FILE__, __LINE__, __func__);
    break;
  case 1:
    
    break;
  case 2:
    // node identification in 2D
    ts = clock();
    create_subdom_graph_2d();
    te = clock();
    tc = tc+te-ts;
    if(mespr == 1) fprintf (out,"Time of creation of subdomain graph %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
    ts = clock();
    create_master_graph();
    te = clock();
    tc = tc+te-ts;
    if(mespr == 1) fprintf (out,"Time of creation of master graph %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
    ts = clock();
    select_minimal_number_2d();
    te = clock();
    tc = tc+te-ts;
    if(mespr == 1) fprintf (out,"Time of selection of minimal number of fixing nodes %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
    if(condfixing == automatic || condfixing == userdef){
      if(myrank == 0){
	ts = clock();
	set_curves();
	te = clock();
	tc = tc+te-ts;
	if(mespr == 1) fprintf (out,"Time of creation of subgraph for adding further fixing nodes %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
	ts = clock();
	if(mespr == 1) print_info_curves();
	te = clock();
	tc = tc+te-ts;
	ts = clock();
       if(condfixing == automatic){
	 ts = clock(); 
	 select_optimum_2d();
	 te = clock();
	 tc = tc+te-ts;
	 if(mespr == 1) fprintf (out,"Time of selection of optimal number of fixing nodes %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
	 // optimum zhruba 5% z hranicnich
       }
       ts = clock(); 
       add_fixings();
       te = clock();
       tc = tc+te-ts;
       if(mespr == 1) fprintf (out,"Time of adding further fixing nodes %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
       reduce_flag = 0;
      }
    }
    else{
      ts = clock(); 
      check_minimal_number_2d();
      te = clock();
      tc = tc+te-ts;
      if(mespr == 1) fprintf (out,"Time of check of minimal number of fixing nodes %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
    }
    ts = clock(); 
    assemble_ltg_with_fixings(ltg);
    te = clock();
    tc = tc+te-ts;
    if(mespr == 1) fprintf (out,"Time of assembling of new ltg array %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
    if(mespr == 1) fprintf (out,"Whole time of fixing node detection  %le s\n",(tc)/(double)CLOCKS_PER_SEC);    
    fflush(out);
    break;
  case 3:
    //node identification in 3d
    switch(condfixing){
    case nocondconer:{
      ts = clock();
      compute_statistics_of_multiplicity();
      te = clock();
      tc = tc+te-ts;
      if(mespr == 1) fprintf (out,"Time of creation of statistics %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
      minSelStrat = 2;
      switch(minSelStrat){
      case 1:
	if(mespr == 1) fprintf(out,"\nFixing nodes will be select by algorithm based on statictics\n");	
	if(mespr == 1) fprintf(stdout,"\nFixing nodes will be select by algorithm based on statictics\n");
	ts = clock();
	select_fixing_nodes_on_master();
	te = clock();
	tc = tc+te-ts;
	if(mespr == 1) fprintf (out,"Time of selection of fixing nodes on master %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
	ts = clock();
	check_minimal_number_3d();
	te = clock();
	tc = tc+te-ts;
	if(mespr == 1) fprintf (out,"Time of check of minimal number of fixing nodes %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
	break;
      case 2:
	if(mespr == 1) fprintf(out,"\nFixing nodes will be select by algortithm based on  topology\n");
	if(mespr == 1) fprintf(stdout,"\nFixing nodes will be select by algortithm based on  topology\n");
	ts = clock();
	create_subdom_graph_3d ();
	te = clock();
	tc = tc+te-ts;
	if(mespr == 1) fprintf (out,"Time of creation of subdomain graph %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
	ts = clock();
	create_master_graph();
	te = clock();
	tc = tc+te-ts;
	if(mespr == 1) fprintf (out,"Time of creation of  master graph %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
	ts = clock();
	select_minimal_number_3d();
	te = clock();
	tc = tc+te-ts;
	if(mespr == 1) fprintf (out,"Time of selection of minimal number of fixing nodes %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
	ts = clock();
	check_minimal_number_3d();
	te = clock();
	tc = tc+te-ts;
	if(mespr == 1) fprintf (out,"Time of check of minimal number of fixing nodes %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
	break;
      }
      ts = clock();
      assemble_ltg_with_fixings(ltg);
      te = clock();
      tc = tc+te-ts;
      if(mespr == 1) fprintf (out,"Time of assembling of new ltg array %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
      if(mespr == 1) fprintf (out,"Whole time of fixing node detection  %le s\n",(tc)/(double)CLOCKS_PER_SEC);    
      fflush(out);
      break;
    }
    case automatic:{
      
      break;
    }
    case userdef:{
      ts = clock();
      create_subdom_graph_3d();
      te = clock();
      tc = tc+te-ts;
      if(mespr == 1) fprintf (out,"Time of creation of subdomain graph %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
      ts = clock();
      create_master_graph();
      te = clock();
      tc = tc+te-ts;
      if(mespr == 1) fprintf (out,"Time of creation of master graph %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
      ts = clock();
      select_minimal_number_3d();
      te = clock();
      tc = tc+te-ts;
      if(mespr == 1) fprintf (out,"Time of selection of minimal number of fixing nodes %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
      ts = clock();
      check_minimal_number_3d();
      te = clock();
      tc = tc+te-ts;
      if(mespr == 1) fprintf (out,"Time of check of minimal number of fixing nodes %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
      if(methodcondcor ==  curvecond || methodcondcor == cursurf) {
	if(myrank == 0){
	  ts = clock();
	  reduce_master_graph_curve();
	  te = clock();
	  tc = tc+te-ts;
	  if(mespr == 1) fprintf (out,"Time of reduction of master graph to curve graph %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
	  ts = clock();
	  set_curves_3D();
	  te = clock();
	  tc = tc+te-ts;
	  if(mespr == 1) fprintf (out,"Time of creation of subgraph for adding further fixing nodes %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
	  ts = clock();
	  if(mespr == 1) print_info_curves();
	  te = clock();
	  tc = tc+te-ts;
	  ts = clock();
	  add_fixings();
	  te = clock();
	  tc = tc+te-ts;
	  if(mespr == 1) fprintf (out,"Time of adding further fixing nodes %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
 	}
      }
      if (methodcondcor ==  surfacecond || methodcondcor == cursurf){
	if(myrank == 0){
	  ts = clock();
	  reduce_master_graph_surface();  
	  if(mespr == 1) print_info_surfaces();
	  te = clock();
	  tc = tc+te-ts;
	  if(mespr == 1) fprintf (out,"Time of reduction of master graph to surface graph %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
	}
	if(typecondsurf == max_tria){
	  ts = clock();
	  add_max_triangle_on_surface();
	  te = clock();
	  tc = tc+te-ts;
	  if(mespr == 1) fprintf (out,"Time of adding further fixing nodes into max triangle %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
	}
	else{
	  ts = clock();
	  select_centers_of_surfaces();
	  te = clock();
	  tc = tc+te-ts;
	  if(mespr == 1) fprintf (out,"Time of selection of surface centres %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
	  if(typecondsurf  == n_mark || typecondsurf == chose_ring || typecondsurf == chose_max_ring){
	    if(myrank == 0){
	      ts = clock();
	      mark_surfnodes();
	      te = clock();
	      tc = tc+te-ts;
	    }
	    // nova fukce
	    if(typecondsurf  == n_mark || typecondsurf == chose_ring){
	      order_selected_ring_nodes();
	    }
	  }
	  if(myrank == 0){
	    ts = clock();
	    if(mespr == 1) fprintf (out,"Time of selection of marking of surface members %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
	    add_fixings();
	    te = clock();
	    tc = tc+te-ts;
	    if(mespr == 1) fprintf (out,"Time of adding further fixing nodes %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    	  
 	  }
 	}
      }
      ts = clock();
      assemble_ltg_with_fixings(ltg);
      te = clock();
      tc = tc+te-ts;
      ts = clock();
      if(mespr == 1) fprintf (out,"Time of assembling of new ltg array %le s\n",(te-ts)/(double)CLOCKS_PER_SEC);    
      if(mespr == 1) fprintf (out,"Whole time of fixing node detection  %le s\n",(tc)/(double)CLOCKS_PER_SEC);    
      fflush(out);
      break;
    }
    }
    fflush(out);
    break;
    
  default:{
    par_print_err(myrank,"Spatial dimension is not detected", __FILE__, __LINE__, __func__);
    break;
  }
  }
}


/**
  Function reduce boundary graph and create boundary curves graph
  It is used for 3D meshes only
  @output curnod
  @output curnadjac
  @output curadjac
  @output curidentif
  JB
**/
void fixnodesel::reduce_master_graph_curve()
{
  long i,j,k,m,l;
  long *correduce;

  reduce_flag = 3;
  
  ncurnod = 0;
  for(i = 0; i < tnbn; i++){
    if(multip[i] > 2){
      ncurnod++;
    }
  }

  curnod = new long[ncurnod];  
  curnadjac = new long[ncurnod];
  curadjac = new long* [ncurnod];
  correduce = new long[tnbn]; 
  curidentif = new long[ncurnod];
  k = 0;
  l = 0;
  for(i = 0; i < tnbn; i++){
    if(multip[i] > 2){
      curnadjac[k] = 0;
      curnod[k] = i;
      correduce[i] = k;
      curidentif[k] = coarseidentif[i];
      for(j = 0; j < coarsenadjac[i]; j++){
	if(multip[coarseadjac[i][j]] > 2 ){
	  curnadjac[k]++;
	}
      }
      curadjac[k] = new long[curnadjac[k]];
      m = 0;
      for(j = 0; j < coarsenadjac[i]; j++){
	if(multip[coarseadjac[i][j]] > 2 ){
	  curadjac[k][m] = coarseadjac[i][j];
	  m++;
	}
      }
      k++;
    }
    else{
      correduce[i] = -1;
    }
  }
  //fprintf(out,"ncurnod je %ld\n",ncurnod);    
  
  // rewrite to reduced numbering
  for(i = 0; i < ncurnod; i++){
    for(j = 0; j < curnadjac[i]; j++){
      m = curadjac[i][j];
      if(correduce[m] != -1){
	curadjac[i][j] = correduce[m];
      }
    }
  }
  
  // detection of cross points 
  for( i = 0; i < ncurnod; i++){
    if(curnadjac[i] > 2 && curidentif[i] != 3){
      curidentif[i] = 3;
    }
    if(curnadjac[i] == 1 && curidentif[i] != 3){
      curidentif[i] = 3;
    }
  }
  
  
  //   for(i = 0; i < ncurnod; i++){
  //     fprintf(out,"uzel %ld identif %ld ma %ld sousedu:",i,curidentif[i],curnadjac[i]);
  //     for(j = 0; j < curnadjac[i]; j++){
  //       fprintf(out,"   %ld",curadjac[i][j]);
  //     }
  //     fprintf(out,"\n");
  //   }
  delete []correduce;
  
  
}

/**
  Function reduce boundary graph and create boundary surfaces graph
  It is used for 3D meshes only

  JB
**/

void fixnodesel::reduce_master_graph_surface()
{
  long i,j,k,m,l,n;
  long a,b;
  long auxnsurf;
  long *aux;
  long **auxsurfmembers;
  long *auxnsurfmembers;
  long **auxsurfnodmark;
  
  reduce_flag = 2;

  
  // surface playing groud
  subdomsurf = new long*[nproc-1];
  for(i = 0; i < nproc-1;i++){
    subdomsurf[i] = new long[nproc-i-1];
    for(j = 0; j < nproc-i-1; j++){
      subdomsurf[i][j] = 0;
    }
  }
  
  for(i = 0; i < tnbn; i++){
    if(multip[i] == 2){
      j = glinkdom[i][0];
      k = glinkdom[i][1];
      if(j < k){
	m = k-(j+1);
	subdomsurf[j][m]++;
      }
      else{
	m = j-(k+1);
	subdomsurf[k][m]++;
      }
    }
  }
  
  // for(i = 0; i < nproc-1;i++){
  //       for(j = 0; j < nproc-i-1; j++){
  //         fprintf(out,"%ld %ld subdomsurf %ld\n",i,j,subdomsurf[i][j]);
  //       }
  //     }
  
  auxnsurf = 0;
  for(i = 0; i < nproc-1; i++){
    for(j = 0; j < nproc-i-1; j++){
      if(subdomsurf[i][j] != 0){
	auxnsurf++;
      }
    }
  }
  //fprintf(out,"Number of auxsurfaces is %ld\n",auxnsurf);

  auxsurfmembers = new long *[auxnsurf];
  auxnsurfmembers = new long [auxnsurf];
  m = 0;
  for(i = 0; i < nproc-1; i++){
    for(j = 0; j < nproc-i-1; j++){
      if(subdomsurf[i][j] != 0){
 	auxsurfmembers[m] = new long[subdomsurf[i][j]];
	auxnsurfmembers[m] = 0;
	subdomsurf[i][j] = m;
	m++;
      }
    }
  }
  
  surfnodpoint = new long[tnbn];
  surfnod = new long[tnbn];
  n = 0;
  for(i = 0; i < tnbn; i++){
    if(multip[i] == 2){
      j = glinkdom[i][0];
      k = glinkdom[i][1];
      m = k-(j+1);
      l = subdomsurf[j][m];
      auxsurfmembers[l][auxnsurfmembers[l]]=i;
      surfnodpoint[i]=auxnsurfmembers[l];
      surfnod[i]=l;
      auxnsurfmembers[l]++;
      n++;
    }
    else{
      surfnodpoint[i]= -1;
      surfnod[i] = -1;
    }
  }
  
  for(i = 0; i < nproc-1; i++){
    delete []subdomsurf[i];
  }
  delete []subdomsurf;
  
  // fprintf(out,"Surfnod\n");
  //   for(i = 0; i < tnbn; i++){
  //     if(multip[i] == 2){
  //       fprintf(out,"%ld surfnod %ld surfnodpoint %ld\n",i,surfnod[i],surfnodpoint[i]);
  //     }
  //   }
  
  // for(i = 0; i < auxnsurf; i++){
  //     fprintf(out,"Surface number %ld",i);
  //     fprintf(out,"has %ld nodes on surface\n",auxnsurfmembers[i]);
  //     for(j = 0; j < auxnsurfmembers[i]; j++){
  //       fprintf(out,"    %ld",auxsurfmembers[i][j]);
  //     }
  //     fprintf(out,"\n");
  //   }
  
  // kontrola na sovislost ploch
  l = -2;
  aux = new long[auxnsurf];
  auxsurfnodmark = new long*[auxnsurf];
  for(i = 0; i < auxnsurf; i++){
    //fprintf(out,"plocha %ld\n",i);
    auxsurfnodmark[i] = new long[auxnsurfmembers[i]];
    for(j = 1; j < auxnsurfmembers[i]; j++){
      auxsurfnodmark[i][j] = -1;
    }
    auxsurfnodmark[i][0] = 0;
    n = 1;
    while(n != auxnsurfmembers[i]){
      for(m = 0; m < auxnsurfmembers[i]; m++){
	for(j = 0;j < auxnsurfmembers[i]; j++){
	  if(auxsurfnodmark[i][j] == m){
	    a = auxsurfmembers[i][j];
	    for(k = 0; k < coarsenadjac[a]; k++){
	      b = coarseadjac[a][k];
	      if(surfnod[b] == i && auxsurfnodmark[i][surfnodpoint[b]] == -1){
		auxsurfnodmark[i][surfnodpoint[b]] = m+1;
		n++;
	      }
	    }
	  }
	  if(n == auxnsurfmembers[i]){
	    break;
	  }
	}
	if(n == auxnsurfmembers[i]){
	  break;
	}
      }
      for(j = 0;j < auxnsurfmembers[i]; j++){
	if(auxsurfnodmark[i][j] >= 0){
	  auxsurfnodmark[i][j] = l;
	}
      }
      aux[i] = -1*l-2;
      l--;
      if(n != auxnsurfmembers[i]){
	for(j = 0;j < auxnsurfmembers[i]; j++){
	  if(auxsurfnodmark[i][j] == -1){
	    if(n != auxnsurfmembers[i]-1){
	      auxsurfnodmark[i][j] = 0;
	      n++;
	    }
	    else{
	      auxsurfnodmark[i][j] = l;
	      n++;
	      aux[i] = -1*l-2;
	      l--;
	    }
	    break;
	  }
	}
      }
    }
 
    // fprintf(out,"%ld %ld\n",n,auxnsurfmembers[i]);
    //     for(j = 0;j < auxnsurfmembers[i]; j++){
    //       fprintf(out,"%ld auxsurfnodmark %ld\n",auxsurfmembers[i][j],auxsurfnodmark[i][j]);
    //     }
  }
  
  nsurf = -1*l-2;
  //fprintf(out,"l je %ld\n",nsurf);
  nsurfmembers = new long [nsurf];
  surfmembers = new long* [nsurf];
  for(i = 0; i < nsurf; i++){
    nsurfmembers[i] = 0;
  } 
  if(aux[0] > 1){
    for(j = 0; j < auxnsurfmembers[0]; j++){
      k = -1*auxsurfnodmark[0][j]-2;
      nsurfmembers[k]++;
    }
    for(j = 0; j < aux[i]+1; j++){
      surfmembers[j] = new long[nsurfmembers[j]];
      n = 0;
      for(k = 0; k < auxnsurfmembers[0]; k++){
	a = -1*auxsurfnodmark[0][k]-2;
	if(a == j){
	  b = auxsurfmembers[0][k];
	  surfnod[b] = j;
	  surfnodpoint[b] = n;
	  surfmembers[j][n] = auxsurfmembers[i][k];
	  n++;
	}
      }
    }
  }
  else{
    nsurfmembers[aux[0]] = auxnsurfmembers[0];
    surfmembers[aux[0]] = new long[auxnsurfmembers[0]];
    for(k = 0; k < auxnsurfmembers[0]; k++){
      b = auxsurfmembers[0][k];
      surfnod[b] = aux[0];
      surfnodpoint[b] = k;
      surfmembers[aux[0]][k] = auxsurfmembers[0][k];
    }
  }
  delete []auxsurfnodmark[0];
  
  for(i = 1; i < auxnsurf; i++){
    if(aux[i]-aux[i-1] > 1 ){
      for(j = 0; j < auxnsurfmembers[i]; j++){
	k = -1*auxsurfnodmark[i][j]-2;
	nsurfmembers[k]++;
      }
      for(j = aux[i-1]+1; j < aux[i]+1; j++){
	surfmembers[j] = new long[nsurfmembers[j]];
 	n = 0;
 	for(k = 0; k < auxnsurfmembers[i]; k++){
	  a = -1*auxsurfnodmark[i][k]-2;
 	  if(a == j){
 	    b = auxsurfmembers[i][k];
 	    surfnod[b] = j;
 	    surfnodpoint[b] = n;
	    surfmembers[j][n] = auxsurfmembers[i][k];
 	    n++;
 	  }
 	}
      }
    }
    else{
      nsurfmembers[aux[i]] = auxnsurfmembers[i];
      surfmembers[aux[i]] = new long[nsurfmembers[aux[i]]];
      for(k = 0; k < auxnsurfmembers[i]; k++){
	b = auxsurfmembers[i][k];
	surfnod[b] = aux[i];
	surfnodpoint[b] = k;
	surfmembers[aux[i]][k] = auxsurfmembers[i][k];
      }
    }
    delete []auxsurfnodmark[i];
  }

  surfdom = new long*[nsurf];
  for(i = 0; i < nsurf; i++){
    surfdom[i] = new long[2];
    a = surfmembers[i][0];
    surfdom[i][0]=glinkdom[a][0];
    surfdom[i][1]=glinkdom[a][1];
  }
  
  
  for(i = 0; i < auxnsurf; i++){
    delete []auxsurfmembers[i];
  }
  delete []auxsurfmembers;
  delete []auxnsurfmembers;
  delete []auxsurfnodmark;

    
 delete []aux;
}

/**
   Function establish interface subgraphs in 2D.
   Intreface subgraphs are defined by basic fixing nodes
   obtained from @function select_fixing_nodes_2D
   JB
 **/

void fixnodesel::set_curves( )
{
  long i,j,k;
  long cor,comb;
  long prev,cur,stop,n,m,a,b,z,v;
  long *startnod,*endnod,**corneib,*corpoint;
  long **auxmemb;
  long *nauxmemb;
  
  // compution of number of start nodes and check of array coarseidentif
  cor = 0;
  if(coarseidentif != NULL){
    for( i = 0; i < tnbn; i++){
      if(coarseidentif[i] == 3 ){
	cor++;
      }
    }
    //fprintf(out,"cor je %ld\n",cor);    
    if(cor == 0){
      for(i = 0; i < tnbn; i++){
	if(coarsenadjac[i] == 1 && coarseidentif[i] != 3 ){
	  coarseidentif[i] = 3;
	  cor++;
	}
	if(coarsenadjac[i] > 2 && coarseidentif[i] != 3){
	  coarseidentif[i] = 3;
	  cor++;
	}
      }
    } 
  }
  else{
    coarseidentif = new long[tnbn];
    for(i = 0; i < tnbn; i++){
      if(coarsenadjac[i] == 1){
	coarseidentif[i] = 3;
      }
      if(coarsenadjac[i] > 2){
	coarseidentif[i] = 3;
	}
    }
  }
  //fprintf(out,"cor je %ld\n",cor);    
  // allocation of auxiliarry array and estimation of number of curves
  corneib = new long*[cor];
  corpoint = new long[cor];
  k = 0;
  comb = 0;
  for(i = 0; i < tnbn; i++){
    if(coarseidentif[i] == 3 && coarsenadjac[i] != 0){
      comb += coarsenadjac[i];
      corpoint[k] = i;
      corneib[k] = new long[coarsenadjac[i]];
      for(j = 0; j < coarsenadjac[i]; j++){
	corneib[k][j] = -1;
      }
      k++;
    }
  }
  if(mespr == 1) fprintf(out,"Estimate number of curves %ld\n",comb);

  //for(i = 0; i < cor; i++){
  //fprintf(out,"%ld corpoint %ld\n",i,corpoint[i]);
  //for(j = 0; j < coarsenadjac[corpoint[i]]; j++){
  //fprintf(out,"   %ld ",corneib[i][j]);
  //}
  //fprintf(out,"\n");
  //}
  
  // old estimation of number of curves
  //comb = compute_number_of_combination(cor,2);
  
  
  startnod = new long[comb];
  endnod = new long[comb];
  auxmemb = new long*[comb];
  nauxmemb = new long[comb];
  
  
  m = 0;
  for(i = 0; i < cor; i++){
    cur = corpoint[i];
    for(j = 0; j < coarsenadjac[corpoint[i]]; j++){
      prev = corpoint[i];
      if(corneib[i][j] == -1){
	startnod[m] = corpoint[i];
	//fprintf(out,"startnode %ld  = %ld\n",i,startnod[m]);
	stop = 0;
	n = 0;
	corneib[i][j] = 0;
	cur = coarseadjac[corpoint[i]][j];
	nauxmemb[m] = 1;
	auxmemb[m] = new long[tnbn];
	auxmemb[m][n] = cur;
	n++;
	while(stop != 1){
	  a = coarsenadjac [cur];
	  for(k = 0; k < a; k++){
	    b = coarseadjac[cur][k];
	    if(b != prev){
	      if(coarseidentif[b] == 2){
		auxmemb[m][n]=b;
		n++;
		prev = cur;
		cur = b;
		stop = 0;
	      }
	      if(coarseidentif[b] == 3){
		stop = 1;
		endnod[m] = b;
		nauxmemb[m]= n;
		//fprintf(out,"konecny bod je %ld\n",b);
		// nastaveni stopky
		for(z = 0; z < cor; z++){
		  if(b == corpoint[z]){
		    //fprintf(out,"corpoint je %ld\n",z);
		    for(v = 0; v < coarsenadjac[corpoint[z]]; v++){
		      if(cur == coarseadjac[corpoint[z]][v]){
			corneib[z][v] = 0;
			//fprintf(out,"nastavena nula na %ld %ld\n",cur,v);
			break;
		      }
		    }
		    break;
		  }
		}
	      }
	    }
	  }
	}
	m++;
      }
    }
  }
  //if(mespr == 1) fprintf(out,"Number of curves %ld\n",m);
  
  ncurves = m;
  start = new long[ncurves];
  // list of end nodes of curves
  end = new long[ncurves];
  // list of members of curves
  members = new long*[ncurves];
  // 
  nmembers = new long[ncurves];
  //
  nedges = new long[ncurves];
  for(i = 0; i < ncurves; i++){
    start[i] = startnod[i];
    end[i] = endnod[i];
    nmembers[i] = nauxmemb[i];
    members[i] = new long[nmembers[i]];
    for(j = 0; j < nmembers[i]; j++){
      members[i][j] = auxmemb[i][j];
    }
    nedges[i] = nmembers[i] + 1; 
  }
  
  delete []startnod;
  delete []endnod;
  delete []nauxmemb;
  for(i = 0; i < ncurves; i++){
    delete []auxmemb[i];
  }
  delete []auxmemb;
  delete []corpoint;
  for(i = 0; i < cor; i++){
    delete []corneib[i];
  }
  delete []corneib;
  
  
}

/**
   Function establish interface curve subgraphs in 3D.
   Intreface curve subgraphs are defined by basic fixing nodes
   obtained from @function select_fixing_nodes_3D
JB
**/

void fixnodesel::set_curves_3D( )
{
  long i,j,k;
  long cor,comb;
  long prev,cur,stop,n,m,a,b,z,v;
  long *startnod,*endnod,**corneib,*corpoint;
  long **auxmemb;
  long *nauxmemb;
  
  cor = 0;
  if(curidentif != NULL){
    for( i = 0; i < ncurnod; i++){
      if(curidentif[i] == 3 ){
	cor++;
      }
    }
    //fprintf(out,"cor je %ld\n",cor);    
    if(cor == 0){
      for(i = 0; i < ncurnod; i++){
	if(curnadjac[i] == 1){
	  curidentif[i] = 3;
	}
       }
    } 
  }
  else{
    curidentif = new long[ncurnod];
    for(i = 0; i < ncurnod; i++){
      if(curnadjac[i] == 1){
	curidentif[i] = 3;
      }
    }
  }

  if(cor == 0){
    for(i = 0; i < ncurnod; i++){
      if(curnadjac[i] == 1){
	curidentif[i] = 3;
      }
    }
  }
  
  cor = 0;
  for( i = 0; i < ncurnod; i++){
    if(curidentif[i] == 3 ){
      cor++;
    }
  }
  
  //fprintf(out,"cor je %ld\n",cor);    

  
  corneib = new long*[cor];
  corpoint = new long[cor];
  k = 0;
  comb = 0;
  for(i = 0; i < ncurnod; i++){
    if(curidentif[i] == 3){
      comb += curnadjac[i];
      corpoint[k] = i;
      corneib[k] = new long[curnadjac[i]];
      for(j = 0; j < curnadjac[i]; j++){
	corneib[k][j] = -1;
      }
      k++;
    }
  }
  
  //for(i = 0; i < redcor; i++){
  //fprintf(out,"%ld corpoint %ld\n",i,corpoint[i]);
  //for(j = 0; j < nadjacnod[corpoint[i]]; j++){
  //fprintf(out,"   %ld ",corneib[i][j]);
  //}
  //fprintf(out,"\n");
  //}
  //old estimation of number of curves
  //comb = compute_number_of_combination(cor,2);
  if(mespr == 1) fprintf(out,"Estimate number of curves %ld\n",comb);
   
  startnod = new long[comb];
  endnod = new long[comb];
  auxmemb = new long*[comb];
  nauxmemb = new long[comb];
  
   
  m = 0;
  for(i = 0; i < cor; i++){
    cur = corpoint[i];
    for(j = 0; j < curnadjac[corpoint[i]]; j++){
      prev = corpoint[i];
      if(corneib[i][j] == -1){
	startnod[m] = corpoint[i];
	//fprintf(out,"startovaci bod je %ld\n",startnod[m]);
	stop = 0;
	n = 0;
	corneib[i][j] = 0;
	cur = curadjac[corpoint[i]][j];
	//fprintf(out,"cur = %ld\n",cur);
	if(curidentif[cur] != 3){
	  nauxmemb[m] = 1;
	  auxmemb[m] = new long[ncurnod];
	  auxmemb[m][n] = cur;
	  n++;
	  while(stop != 1 ){
	    a = curnadjac [cur];
	    for(k = 0; k < a; k++){
	      b = curadjac[cur][k];
	      if(b != prev){
		//fprintf(out,"b=%ld\n",b);
		if(curidentif[b] == 2){
		  auxmemb[m][n]=b;
		  n++;
		  prev = cur;
		  cur = b;
		  // fprintf(out,"cur = %ld\n",cur);
		  stop = 0;
		  //break;
		}
		if(curidentif[b] == 3){
		  stop = 1;
		  endnod[m] = b;
		  nauxmemb[m]= n;
		  //fprintf(out,"konecny bod je %ld\n",b);
		  // nastaveni stopky
		  for(z = 0; z < cor; z++){
		    if(b == corpoint[z]){
		      //fprintf(out,"corpoint je %ld\n",z);
		      for(v = 0; v < curnadjac[corpoint[z]]; v++){
			if(cur == curadjac[corpoint[z]][v]){
			  corneib[z][v] = 0;
			  //fprintf(out,"nastavena nula na %ld %ld\n",cur,v);
			  break;
			}
		      }
		      break;
		    }
		  }
		  break;
		}
	      }
	    }
	  }
	  m++;
	}
	else{
	  //fprintf(out,"konecny bod je %ld\n",cur);
	  nauxmemb[m] = 0;
	  endnod[m] = cur;
	  // identfikace cile
	  for(z = 0; z < cor; z++){
	    if(cur == corpoint[z]){
	      //fprintf(out,"corpoint je %ld\n",z);
	      for(v = 0; v < curnadjac[corpoint[z]]; v++){
		if(prev == curadjac[corpoint[z]][v]){
		  corneib[z][v] = 0;
		  //fprintf(out,"nastavena nula na %ld %ld\n",cur,v);
		  break;
		}
	      }
	      break;
	    }
	  }
	}
      }
    }
  }
  //fprintf(out,"Number of curves %ld\n",m);
  
  ncurves = m;
  start = new long[ncurves];
  // list of end nodes of curves
  end = new long[ncurves];
  // list of members of curves
  members = new long*[ncurves];
  // 
  nmembers = new long[ncurves];
  //
  nedges = new long[ncurves];
  for(i = 0; i < ncurves; i++){
    start[i] = startnod[i];
    end[i] = endnod[i];
    nmembers[i] = nauxmemb[i];
    members[i] = new long[nmembers[i]];
    for(j = 0; j < nmembers[i]; j++){
      // prepis z lokalniho krivkoveho cislovani na globalni
      k = curnod[auxmemb[i][j]];
      members[i][j] = k;
    }
    nedges[i] = nmembers[i] + 1; 
  }
  
  delete []startnod;
  delete []endnod;
  delete []nauxmemb;
  for(i = 0; i < ncurves; i++){
    delete []auxmemb[i];
  }
  delete []auxmemb;
  delete []corpoint;
  for(i = 0; i < cor; i++){
    delete []corneib[i];
  }
  delete []corneib;
  
}


/**
   Function prints information about interface curves (interface subgraphs
   in 2D or interface curve subgraphs in 3D)
   JB
 **/
void fixnodesel::print_info_curves()
{
  long i,j;
  
  fprintf(out,"Number of curves %ld\n",ncurves);
  for(i = 0; i < ncurves; i++){
    fprintf(out,"Curve %ld:\n",i+1);
    fprintf(out,"Start %ld\n",start[i]);
    fprintf(out,"End %ld\n",end[i]);
    fprintf(out,"Number of members %ld\n",nmembers[i]);
    fprintf(out,"Members: ");
    for(j = 0; j < nmembers[i]; j++){
      fprintf(out,"   %ld",members[i][j]);
    }
    fprintf(out,"\n");
    fprintf(out,"Number of edges %ld\n",nedges[i]);
  }
}


/**
   Function prints information about interface surfaces 
   (interface surface subgraphs)
   JB
 **/
void fixnodesel::print_info_surfaces()
{
  long i,j;
  fprintf(out,"Number of surfaces is %ld\n",nsurf);
  for(i = 0; i < nsurf; i++){
    fprintf(out,"Surface number %ld lie on subdomain %ld and %ld\n",i,surfdom[i][0],surfdom[i][1]);
    fprintf(out,"and has %ld nodes\nMembers:  ",nsurfmembers[i]);
    for(j = 0; j < nsurfmembers[i]; j++){
      fprintf(out,"    %ld",surfmembers[i][j]);
    }
    fprintf(out,"\n");
  }
}


/**
   Function adds fixing nodes in accroding to user requirement
   JB
 **/
void fixnodesel::add_fixings()
{
  //fprintf(out,"methodcondcor %ld\n",methodcondcor);
  if(methodcondcor ==  curvecond || methodcondcor == cursurf){
    switch(typecondcur){
    case nth_memb:{
      //add_each_nth_member();
      add_nth_member();
      break;
    }
    case all_memb:{
      add_all_members();
      break;
    }
    case rand_memb:{
      add_rand();
      break;
    }
    case centroid_fix:{
      add_centroid();
      break;
    }
    case n_part_curve:{
      add_n_part_curve();
      break;
    }
    case userposdef:{
      add_user_pos_def();
      break;
    }
    }
  }
  if (methodcondcor ==  surfacecond || methodcondcor == cursurf){
    //fprintf(out,"type %ld\n",typecondsurf);
    switch(typecondsurf){
    case rand_memb:{
      add_n_rand_nodes_on_surf();
      break;
    }
    case all_memb:{
      add_all_surfnod();
      break;
    }
    case centroid_fix:{
      add_centroid_surface();
      break;
    }
    case n_mark:{
      add_n_th_mark();
      break;
    }
    case chose_ring:{
      add_rings();
      break;
    }
    case chose_max_ring:{
	add_max_ring();
	break;
    }
    case userposdef:{
      add_user_pos_def();
      break;
    }  
    }
  }
}


/**
   Function randomly adds fixng nodes into curve graphs
   JB
 **/
void fixnodesel::add_rand()
{
  long i,j,k;
  long member;
  long *sel_rand;
  
  // inicialization of pseudorandom numbers
  srand((unsigned int) time(NULL));
  sel_rand = new long[nmembercur];
  for(i = 0; i < ncurves; i++){
    if(mespr == 1) fprintf(out,"On %ld curve will be selected:\n",i);
    if(nmembers[i] != 0){
      sel_rand[0] = rand()%nmembers[i];
      for(j = 1; j < nmembercur; j++){
	member = rand()%nmembers[i];
	// stop = 0;
// 	while(stop == 0){
// 	  m = 0;
// 	  for(k = 0; k < j; k++){
// 	    if(member == sel_rand[k]){
// 	      m++;
// 	    }
// 	  }
// 	  if(m == 0){
// 	    stop = 1;
// 	  }
// 	  else{
// 	    member = rand()%nmembers[i];
// 	    stop = 0;
// 	  }
// 	}
	sel_rand[j] = member;
      }
      for(j = 0; j < nmembercur; j++){
	k = members[i][sel_rand[j]];
	coarseidentif[k] = 3;
	if(mespr == 1) fprintf(out,"   %ld",members[i][sel_rand[j]]);
      }
    }
    if(mespr == 1) fprintf(out,"\n");
  }
  
 
  delete []sel_rand;
  fflush(out);

}

/**
   Function adds each n-th vertex of subgraph as fixing node
   JB
 **/

void fixnodesel::add_each_nth_member()
{
  long i,j,stop,e;
  double zbytek;
  long cycles,member,lnmember;
  
  for(i = 0; i < ncurves; i++){
    if(mespr == 1) fprintf(out,"On curve %ld will be selected:",i);
    lnmember = nmembercur;
    if(nmembers[i] != 0){
      if(lnmember > nedges[i]){
	//fprintf(out,"nmember is too big\n");
	lnmember = nedges[i]-1;
      }
      zbytek = nedges[i] % lnmember;
      if(zbytek == 0){
	cycles = nedges[i] / lnmember;
      }
      else{
	// nova volba nmember
	stop = 0;
	e = 0;
	while (stop == 0){
	  zbytek = nedges[i] % lnmember;
	  if( zbytek != 0){
	    //fprintf(out,"zbytek = %lf\n",zbytek);
	    if(e == 0){
	      lnmember++;
	    }
	    if(e == 1){
	      lnmember--;
	    }
	    //fprintf(out,"lnmember = %ld\n",lnmember);
	    stop = 0;
	  }
	  else{
	    cycles = nedges[i] / lnmember;
	    //fprintf(out,"nova volba nmember = %ld\n",lnmember);
	    stop = 1;
	  }
	  if( nedges[i] == lnmember){
	    e = 1;
	    lnmember = nmembercur;
	    stop = 0;
	  }
	  if(lnmember == 1){
	    stop = 1;
	    e = 2;
	  }
	}
      }
      if(e == 2){
	//fprintf(out,"%ld je prvocislo, nelze delit - nelze vybrat kazdy n-ty\nbude vybran n-ty clen\n",nedges[i]);
	member = 0;
	for(j = 0; j < nmembers[i]; j++){
	  member++;
	  if(member == nmembercur){
	    //  fprintf(out,"vybran bude uzel cislo %ld\n",members[i][j]);
	    member = 0;
	    coarseidentif[members[i][j]] = 3;
	  }
	}
      }
      //fprintf(out,"cycles = %ld\n",cycles);
      for(j = 1; j < cycles; j++){
	member = j*lnmember-1;
	if(mespr == 1) fprintf(out,"   %ld",members[i][member]);
	coarseidentif[members[i][member]] = 3;
      }
    }
    if(mespr == 1) fprintf(out,"\n");
  }
  fflush(out);
}


/**
   Function adds each n-th vertex of subgraph as fixing node
   JB
 **/
void fixnodesel::add_nth_member()
{
  long i,j;
  long member;
  //fprintf(out,"nmembercur %ld\n",nmembercur);
  for(i = 0; i < ncurves; i++){
    if(condfixing == automatic){
      nmembercur = automember[i];
    }
    
    if(mespr == 1) fprintf(out,"On curve %ld will be selected:",i);
    if(nmembers[i] != 0){
      if(nmembercur >= nmembers[i]){
	//fprintf(out,"nmember is too big\n");
      }
      else{
	member = 0;
	for(j = 0; j < nmembers[i]; j++){
	  member++;
	  if(member == nmembercur){
	    if(mespr == 1) fprintf(out,"   %ld",members[i][j]);
	    member = 0;
	    coarseidentif[members[i][j]] = 3;
	  }
	}
      }
    }
    if(mespr == 1) fprintf(out,"\n");
  }
  fflush(out);
}


/**
   Function adds all vertices of subgraph as fixing nodes
   JB
 **/

void fixnodesel::add_all_members( )
{
  long i,j;
  for(i = 0; i < ncurves; i++){
    if(mespr == 1) fprintf(out,"On curve %ld will be selected:",i);
    for(j = 0; j < nmembers[i]; j++){
      if(mespr == 1) fprintf(out," %ld",members[i][j]);
      coarseidentif[members[i][j]] = 3;
    }
    if(mespr == 1) fprintf(out,"\n");
  }
}


/**
   Function adds centre of subgraph as fixing node
   JB
 **/

void fixnodesel:: add_centroid()
{
  long i,k;
  long select;
  for(i = 0; i < ncurves; i++){
    if((nedges[i]%2) == 0){
      select = nedges[i]/2;
      select--;
    }
    else{
      k = nedges[i]+1;
      select = k/2;
      select--;
    }
    if(mespr == 1) fprintf(out,"On curve %ld will be selected %ld\n",i,members[i][select]);
    coarseidentif[members[i][select]] = 3;
  }

}


// /**
//    Function adds integral points of subgraphs as fixing node
//    JB
//  **/

// void fixnodesel:: add_intpoint()
// {
//   // numbers of nodes which will be added
//   long i;
//   double dpos;
//   long lpos;
//   // intpointo - number of integration points
  
//   switch(intpointocur){
//   case 1:
//     for(i = 0; i < ncurves; i++){
//       fprintf(out,"On curve %ld will be selected:",i);
//       // 1st
//       dpos = 0.5*nedges[i];
//       round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       fprintf(out,"\n");
//     }
//     break;
//   case 2:
//     for(i = 0; i < ncurves; i++){
//       fprintf(out,"On curve %ld will be selected:",i);
//       // 1st
//       dpos = 0.2*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       //fprintf(out,"Curve %ld:vybran bude uzel cislo %ld\n",i,members[i][lpos]);
//       // 2nd
//       dpos = 0.8*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       fprintf(out,"\n");
//     }
//     break;
//   case 3:
//     for(i = 0; i < ncurves; i++){
//       fprintf(out,"On curve %ld will be selected:",i);
//       // 1st
//       dpos = 0.1*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 2nd
//       dpos = 0.5*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 3rd
//       dpos = 0.9*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       fprintf(out,"\n");
//     }
//     break;
//   case 4:
//     for(i = 0; i < ncurves; i++){
//       fprintf(out,"On curve %ld will be selected:",i);
//       // 1st
//       dpos = 0.07*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 2nd
//       dpos = 0.3*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 3rd
//       dpos = 0.7*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 4th
//       dpos = 0.93*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       fprintf(out,"\n");
//     }
//     break;
//   case 5:
//     for(i = 0; i < ncurves; i++){
//       fprintf(out,"On curve %ld will be selected:",i);
//       // 1st
//       dpos = 0.07*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 2nd
//       dpos = 0.23*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 3rd
//       dpos = 0.5*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 4th
//       dpos = 0.77*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 5th
//       dpos = 0.93*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       fprintf(out,"\n");
//     }
//     break;
//   case 6:
//     for(i = 0; i < ncurves; i++){
//       fprintf(out,"On curve %ld will be selected:",i);
//       //1st
//       dpos = 0.03*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 2nd
//       dpos = 0.17*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 3rd
//       dpos = 0.38*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 4rd
//       dpos = 0.62*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 5th
//       dpos = 0.83*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 6th
//       dpos = 0.97*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       fprintf(out,"\n");
//     }
//     break;
//   case 7:
//     fprintf(out,"case 7\n");
//     for(i = 0; i < ncurves; i++){
//       fprintf(out,"On curve %ld will be selected:",i);
//       //1st
//       dpos = 0.03*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 2nd
//       dpos = 0.13*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 3rd
//       dpos = 0.30*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 4th
//       dpos = 0.50*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 5th
//       dpos = 0.70*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 6th
//       dpos = 0.87*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 7th
//       dpos = 0.97*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       fprintf(out,"\n");
//     }
//     break;
//   case 8:
//     for(i = 0; i < ncurves; i++){
//       fprintf(out,"On curve %ld will be selected:",i);
//       //1st
//       dpos = 0.02*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 2nd
//       dpos = 0.10*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 3rd
//       dpos = 0.24*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 4th
//       dpos = 0.41*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 5th
//       dpos = 0.59*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 6th
//       dpos = 0.76*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 7th
//       dpos = 0.90*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       // 8th
//       dpos = 0.98*nedges[i];
//       dpos = round_it(dpos);
//       lpos = (long) dpos;
//       lpos--;
//       coarseidentif[members[i][lpos]] = 3;
//       fprintf(out,"   %ld",members[i][lpos]);
//       fprintf(out,"\n");
//     }
//     break;
//   }
  
//   fflush(out);
// }


/**
   Function adds terminuses of subpaths of subgraphs as fixing nodes
   JB
 **/

void fixnodesel::add_n_part_curve()
{
  long i,j,k;
  long smember,stop,e,lnmember,lnedges;
  double zbytek;
  
  for(i = 0; i < ncurves; i++){
    if(condfixing == automatic){
      nmembercur = automember[i];
    }
    if(mespr == 1) fprintf(out,"On curve %ld will be selected:",i);
    if(nedges[i] > nmembercur ){
      if(nedges[i]%nmembercur == 0) {
	smember = nedges[i]/nmembercur;
      } 
      else{
	lnmember = nmembercur;
	lnedges = nedges[i];
	// nova volba nmember
	stop = 0;
	e = 0;
	while (stop == 0){
	  zbytek = lnedges % lnmember;
	  if( zbytek != 0){
	    //fprintf(out,"zbytek = %lf\n",zbytek);
	    if(e == 0){
	      lnmember--;
	    }
	    if(e == 1){
	      lnmember++;
	    }
	    //fprintf(out,"lnmember = %ld\n",lnmember);
	    stop = 0;
	  }
	  else{
	      smember = lnedges / lnmember;
	      //fprintf(out,"nova volba nmember = %ld\n",lnmember);
	      stop = 1;
	  }
	  if( lnedges == lnmember){
	    e = 1;
	    lnmember = nmembercur;
	    stop = 0;
	  }
	  if(lnmember == 1){
	    stop = 1;
	    e = 2;
	    //fprintf(out,"%ld je prvocislo, nelze delit - nelze vybrat kazdy n-ty\nbude upraven pocet hran na nedges+1 %ld\n",lnedges);
	    lnedges++;
	    lnmember = nmembercur;
	  }
	}
      }
      //fprintf(out,"cycle %ld nmembers %ld \n",smember,nmembers[i]);
      if(e != 2){
	k = 0;
	for(j = 0; j < nmembers[i];j++){
	  k++;
	  if(k == smember){
	    k = 0;
	    if(mespr == 1) fprintf(out,"   %ld",members[i][j]);
	    coarseidentif[members[i][j]] = 3;
	  }
	}
      }
    }
    if(mespr == 1) fprintf(out,"\n");
  }
  
}

/**
   Function adds fixing nodes into user defined positions
   JB
 **/
void fixnodesel::add_user_pos_def()
{
  long i,j,k,select,numbering;
  if(methodcondcor ==  curvecond || methodcondcor == cursurf){
    numbering =1;
    for(i = 0; i < ncurves; i++){
      if(mespr == 1) fprintf(out,"On curve %ld will be selected\n",i);
      coarseidentif[start[i]] = 2;
      coarseidentif[end[i]] = 2;
      for(j = 0; j < nmembers[i]; j++){
	coarseidentif[members[i][j]] = 2;
      }
      //fprintf(out,"Selected nodes on %ld curve, all others fixings will be removed\n",i);
      for(j = 0; j < nuserdefnod; j++){
	select = userdefnod[j]-1;
	if(select == -1 || numbering == 0){ 
	  select++;
	  numbering = 0;
	}
	//fprintf(out," %ld",select);
	if(select < nmembers[i]){
	  k = members[i][select];
	  coarseidentif[k] = 3;
	  if(mespr == 1) fprintf(out," %ld",k);
	}
	else{
	  select = nmembers[i]-1;
	  k = members[i][select];
	  coarseidentif[k] = 3;
	  if(mespr == 1) fprintf(out," %ld",k);
	  break;
	}
      }
      if(mespr == 1) fprintf(out,"\n");
    }
  }
  if (methodcondcor ==  surfacecond || methodcondcor == cursurf){
    numbering =1;
    for(i = 0; i < nsurf; i++){
      if(mespr == 1) fprintf(out,"On surface %ld will be selected\n",i);
      for(j = 0; j < nsurfmembers[i]; j++){
	coarseidentif[surfmembers[i][j]] = 2;
      }
      //fprintf(out,"Selected nodes on %ld curve, all others fixings will be removed\n",i);
      for(j = 0; j < nuserdefnod; j++){
	select = userdefnod[j]-1;
	if(select == -1 || numbering == 0){ 
	  select++;
	  numbering = 0;
	}
	//fprintf(out," %ld",select);
	if(select < nsurfmembers[i]){
	  k = surfmembers[i][select];
	  coarseidentif[k] = 3;
	  if(mespr == 1) fprintf(out," %ld",k);
	}
	else{
	  select = nmembers[i]-1;
	  k = surfmembers[i][select];
	  coarseidentif[k] = 3;
	  if(mespr == 1) fprintf(out," %ld",k);
	  break;
	}
      }
      if(mespr == 1) fprintf(out,"\n");
    }
  }
  
}

/**
   Function randomly adds fixing nodes on surfaces
   JB
 **/
void fixnodesel::add_n_rand_nodes_on_surf()
{

  long i,j,k;
  long member;
  long *sel_rand;
  
  // inicialization of pseudorandom numbers
  srand((unsigned int) time(NULL));
  sel_rand = new long[nmembersurf];
  for(i = 0; i < nsurf; i++){
    if(mespr == 1) fprintf(out,"On %ld surface will be selected:\n",i);
    if(nsurfmembers[i] != 0){
      sel_rand[0] = rand()%nsurfmembers[i];
      for(j = 1; j < nmembersurf; j++){
 	member = rand()%nsurfmembers[i];
// 	stop = 0;
// 	while(stop == 0){
// 	  m = 0;
// 	  for(k = 0; k < j; k++){
// 	    if(member == sel_rand[k]){
// 	      m++;
// 	    }
// 	  }
// 	  if(m == 0){
// 	    stop = 1;
// 	  }
// 	  else{
// 	    member = rand()%nsurfmembers[i];
// 	    stop = 0;
// 	  }
// 	}
 	sel_rand[j] = member;
      }
      for(j = 0; j < nmembersurf; j++){
	k = surfmembers[i][sel_rand[j]];
	coarseidentif[k] = 3;
	if(mespr == 1) fprintf(out,"   %ld",k);
      }
    }
    if(mespr == 1) fprintf(out,"\n");
  }
  
 
  delete []sel_rand;
  fflush(out);
}

/**
   Function adds all nodes on surfaces as fixing nodes
   JB
 **/
void fixnodesel::add_all_surfnod()
{
  long i,j;
  for(i = 0; i < nsurf; i++){
    if(mespr == 1) fprintf(out,"On surface %ld will be selected:",i);
    for(j = 0; j < nsurfmembers[i]; j++){
      if(mespr == 1) fprintf(out," %ld",surfmembers[i][j]);
      coarseidentif[surfmembers[i][j]] = 3;
    }
    if(mespr == 1) fprintf(out,"\n");
  }
}

/**
   Function adds fixing nodes into centres of interface surface subgraphs
   The centre of interface surface subgraphs is computed as a average form 
   geometrical coordinates of vertices of such subgraph
   JB
**/
void fixnodesel::add_centroid_surface()
{
  long i,j,a;
  for(i = 0; i < nsurf; i++){
    if(mespr == 1) fprintf(out,"On surface %ld will be selected:",i);
    for(j = 0; j < nsurfcentres[i]; j++){
      a  =  surfcenters[i][j];
      coarseidentif[a] = 3;
      if(mespr == 1) fprintf(out," %ld",a);
    }
    if(mespr == 1) fprintf(out,"\n");
  }

}
/**
   Function adds fixing nodes on n-th ring, which is created
   with vertices with same length of walk form centre to appropriate
   node
   JB
 **/
void fixnodesel::add_n_th_mark()
{
  long i,j,a;
  for(i = 0; i < nsurf; i++){
    if(mespr == 1) fprintf(out,"On surface %ld will be selected:",i);
    if(nsurfmembers[i] != 0){
      for(j = 0; j < nsurfmembers[i]; j++){
	if((surfnodmark[i][j]%nmembersurf) == 0){
	  a = surfmembers[i][j];
	  coarseidentif[a] = 3;
	  if(mespr == 1) fprintf(out," %ld",surfmembers[i][j]);
	}
      }
    }
    if(mespr == 1) fprintf(out,"\n");
  }
}

/**
   odstranit
   JB
 **/
void fixnodesel::add_int_points_surface()
{
  
}

/**

   JB
**/
void fixnodesel::order_selected_ring_nodes()
{
  long i,j,a,k,max;
  long nring;
  long *ringnodes;
  //  long *ringmark;
  long **pointer;
  //  bool stop;
  MPI_Status stat;

  
  ringnodes = new long[nbn];
  // nmebersurf is each n-th ring
  if(myrank == 0){
    pointer = new long*[nproc];
    for(i = 0; i< nproc;  i++){
      pointer[i] = new long[nbn];
      for (j = 0; j < nbn; j++){
	pointer[i][j] = -1;
      }
    }
    long aa,bb;
    for(i = 0; i < nsurf; i++){
      max = 0;
      aa = surfdom[i][0];
      bb = surfdom[i][1];
      for(j = 0; j < nsurfmembers[i]; j++){
	if(max < surfnodmark[i][j]){
	  max = surfnodmark[i][j];
	}
      }
      if((max%nmembersurf) == 0){
	nring = max/nmembersurf;
      }
      a = nmembersurf;
      for (j = 0; j < nring; j++){
	for(k = 0; k < nsurfmembers[i]; k++){
	  if(surfnodmark[i][k] == a){
	    // k uzel na plose i -> m uzel na podoblasti n
	    pointer[aa][surfnodpoint[surfmembers[i][k]]] = i;
	    pointer[bb][surfnodpoint[surfmembers[i][k]]] = i;
	  }
	  a+=nmembersurf;
	}
      }
    }
    // send
    for(i = 1; i < nproc; i++){
      MPI_Send (pointer[i],nbnd[i],MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    for(i = 0; i < nbnd[0]; i++){
      ringnodes[i] = pointer[0][i];
    }
  }
  else{
    // recieve
    MPI_Recv (ringnodes,nbn,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  // zpracovani a odeslani zpet

  long *aux;
  max = 0;
  for(i = 0; i < nbn; i++){
    if(max < ringnodes[i]){
      max = ringnodes[i];
    }
  }

  max++;
  aux = new long[max];
  for(i = 0; i < max; i++){
    aux[i] = 0;
  }
  
  for(i = 0; i < nbn; i++){
    if(ringnodes[i] != -1){
      j = ringnodes[i];
      aux[j]++;
    }
  }
  double *nodelength;
  double *nodeangle;
  nodelength = new double[nbn];
  nodeangle = new double[nbn];
  k = 0;
  for(i = 0; i < max; i++){
    if(aux[i] != 0){
      for(j = 0; j < nbn; j++){
	if(ringnodes[j] == i){
	  nodelength[j] = (top->gnodes[lnbn[j]].x-realcg[k][0])*(top->gnodes[lnbn[j]].x-realcg[k][0])+(top->gnodes[lnbn[j]].y-realcg[k][1])*(top->gnodes[lnbn[j]].y-realcg[k][1])+(top->gnodes[lnbn[j]].z-realcg[k][2])*(top->gnodes[lnbn[j]].z-realcg[k][2]);
	  nodeangle[j] = (top->gnodes[lnbn[j]].x-realcg[k][0])*(top->gnodes[lnbn[j]].x-realcg[k][0])+(top->gnodes[lnbn[j]].y-realcg[k][1])*(top->gnodes[lnbn[j]].y-realcg[k][1])+(top->gnodes[lnbn[j]].z-realcg[k][2])*(top->gnodes[lnbn[j]].z-realcg[k][2]);
	}
      }

      
      k++;
    }
  }
  
  
}

/**
   Function adds fixing nodes into vertices which lay on user defined ring
   JB
**/

void fixnodesel::add_rings()
{
  long i,j,k,a;
  for(i = 0; i < nsurf; i++){
    if(mespr == 1) fprintf(out,"On surface %ld will be selected:",i);
    if(nring != 0){
      for(j = 0; j < nring; j++){
	for(k = 0; k < nsurfmembers[i]; k++){
	  if(surfnodmark[i][k] == ring[j]){
	    a = surfmembers[i][j];
	    coarseidentif[a] = 3;
	    if(mespr == 1) fprintf(out," %ld",surfmembers[i][k]);
	  }
	}
	if(mespr == 1) fprintf(out,"\n");
      }
    }
    if(mespr == 1) fprintf(out,"\n");
  } 

}

/**
   Function adds fixing nodes on ring with maximal mark
 **/
void fixnodesel::add_max_ring()
{
  long i,j,a;
  long max;
  for(i = 0; i < nsurf; i++){
    max = 0;
    for(j = 0; j < nsurfmembers[i]; j++){
      if(max < surfnodmark[i][j]){
	max = surfnodmark[i][j];
      }
    }
    if(mespr == 1) fprintf(out,"On surface %ld will be selected:",i);
    for(j = 0; j < nsurfmembers[i]; j++){
      if(surfnodmark[i][j] == max || surfnodmark[i][j] == max-1){
	a = surfmembers[i][j];
	coarseidentif[a] = 3;
	if(mespr == 1) fprintf(out," %ld",surfmembers[i][j]);
      }
    }
    if(mespr == 1) fprintf(out,"\n");
  }
  
}
/**
   Function adds fixing nodes on vertices of interface surface subgraphs
   which create triagle with maximal area.
   JB
**/

void fixnodesel::add_max_triangle_on_surface()
{
  long i,j,k,l,max;
  long buffsize,*buff,*aux;
  MPI_Status stat;
  // double t1,t2,tc;

  double maxradius,radius,maxarea,area;
  long **auxnodes;
  long *vcomb;
  long subdomnsurf;
  long *subdomnsurfmem;
  long **subdomsurfmem;

  buffsize = maxnbn;
  buff = new long[buffsize];
  if(myrank == 0){
    for(i = 1; i < nproc; i++){
      for(j = 0; j < nbnd[i]; j++){
	buff[j] = surfnod[cnbn[i][j]];
      }
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
      
    }
    for(j = 0; j < nbnd[0]; j++){
      buff[j] = surfnod[cnbn[0][j]];
    }
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);      
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  max = 0;
  for(i = 0; i < nbn; i++){
    if(max < buff[i]){
      max = buff[i];
    }
  }
  //fprintf(out,"max je %ld\n",max);

  max++;
  aux = new long[max];
  for(i = 0; i < max; i++){
    aux[i] = 0;
  }
  
  for(i = 0; i < nbn; i++){
    if(buff[i] != -1){
      j = buff[i];
      aux[j]++;
    }
  }
  
  subdomnsurf = 0;
  for(i = 0; i < max; i++){
    if(aux[i] != 0){
      subdomnsurf++;
    }
  }
  
  subdomnsurfmem  = new long[subdomnsurf];
  subdomsurfmem  = new long*[subdomnsurf];
  for(i = 0; i < subdomnsurf; i++){
    subdomnsurfmem[i] = 0;
  }
  
  j = 0;
  for(i = 0; i < max; i++){
    if(aux[i] != 0){
      subdomsurfmem[j] = new long[aux[i]];
      aux[i] = j;
      j++;
    }
  }
  
  for(i = 0; i < nbn; i++){
    if(buff[i] != -1){
      j = buff[i];
      k = aux[j];
      subdomsurfmem[k][subdomnsurfmem[k]]=i;
      subdomnsurfmem[k]++;
      buff[i] = 0;
    }
  }
  delete[]aux;
  
//   if(mespr == 1){
//     fprintf(out,"Subdomain %d has %ld surfaces\n",myrank,subdomnsurf);
//     for(i = 0; i < subdomnsurf; i++){
//       fprintf(out,"Surface %ld has %ld nodes\n",i,subdomnsurfmem[i]);
//       for(j = 0; j < subdomnsurfmem[i]; j++){
//         fprintf(out,"   %ld",lgnbn[subdomsurfmem[i][j]]);
//       }
//       fprintf(out,"\n");
//     }
//   }
  
  max = 0;
  for(i = 0; i < subdomnsurf; i++){
    if(max < subdomnsurfmem[i]){
      max = subdomnsurfmem[i];
    }
  }
  
  vcomb = new long[3];
  for(i = 0; i < subdomnsurf; i++){
    vcomb[0] = lnbn[subdomsurfmem[i][0]];
    maxradius = 0.0;
    for(j = 0; j < subdomnsurfmem[i]; j++){
      vcomb[1] = lnbn[subdomsurfmem[i][j]];
      radius = compute_length_of_vector(vcomb[0],vcomb[1]);
      //fprintf(out,"radius je %lf\n\n",radius);
      if(radius > maxradius){
	maxradius = radius;
	k = j;
      }
    }
    vcomb[0] = lnbn[subdomsurfmem[i][k]];
    //buff[subdomsurfmem[i][k]]=1;

    maxradius = 0.0;
    for(j = 0; j < subdomnsurfmem[i]; j++){
      vcomb[1] = lnbn[subdomsurfmem[i][j]];
      radius = compute_length_of_vector(vcomb[0],vcomb[1]);
      if(radius > maxradius){
	k = j;
	maxradius = radius;
      }
    }
    vcomb[0] = lnbn[subdomsurfmem[i][k]];
    buff[subdomsurfmem[i][k]]=1;
    
    maxradius = 0.0;
    for(j = 0; j < subdomnsurfmem[i]; j++){
      vcomb[1] = lnbn[subdomsurfmem[i][j]];
      radius = compute_length_of_vector(vcomb[0],vcomb[1]);
      if(radius > maxradius){
	k = j;
	maxradius = radius;
      }
    }
    vcomb[1] = lnbn[subdomsurfmem[i][k]];
    buff[subdomsurfmem[i][k]]=1;
    
    
    maxarea = 0.0;
    for(j = 0; j < subdomnsurfmem[i]; j++){
      if(vcomb[0] !=  lnbn[subdomsurfmem[i][j]] && vcomb[1] !=  lnbn[subdomsurfmem[i][j]]){
	vcomb[2] = lnbn[subdomsurfmem[i][j]];
	area = compute_area_of_triangle(vcomb[0],vcomb[1],vcomb[2]);
	if(maxarea < area){
	  maxarea = area;
	  k = j;
	}
      }
    }
    buff[subdomsurfmem[i][k]]=1;
     
  }
  delete []vcomb;
  
  
  if(myrank == 0){
    auxnodes = new long*[nproc];
    auxnodes[0] = new long[nbnd[0]];
    for(j = 0; j < nbnd[0]; j++){
      auxnodes[0][j] = buff[j];
    }
    for(i = 1; i < nproc; i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k = domproc[stat.MPI_TAG];
      auxnodes[k] = new long[nbnd[k]];
      for(j = 0; j < nbnd[k]; j++){
	auxnodes[k][j] = buff[j];
      }
    }
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete []buff;
  
  if(myrank == 0){
    for(i = 0; i < nproc; i++){
      for(j = 0; j < nbnd[i]; j++){
	if(auxnodes[i][j] == 1){
	  k = cnbn[i][j];
	  l = surfnod[k];
	  if(glinkdom[k][1] <= i){
	    coarseidentif[k] = 3;
	  }
	}
      }
    }
    
    for(i = 0; i < nsurf; i++){
      if(mespr == 1) fprintf(out,"On surface %ld will be selected:",i);
      for(j = 0; j < nsurfmembers[i]; j++){
	k = surfmembers[i][j];
	if(coarseidentif[k] == 3){
	  if(mespr == 1) fprintf(out," %ld",surfmembers[i][j]);
	}
      }
      if(mespr == 1) fprintf(out,"\n");
    }
    for(i = 0; i < nproc; i++){
      delete []auxnodes[i];
    }
    delete []auxnodes;
  }
}



