#include "axisrotrec.h"
#include "probdesc.h"
#include "mechtop.h"
#include "global.h"
#include "globmat.h"
#include "iotools.h"
#include "vector.h"
#include "matrix.h"


axisrotrec::axisrotrec() : l(3)
{
  dax = noaxdef;
  a[0] = a[1] = a[2] = b[0] = b[1] = b[2] = 0.0;
  u1[0] = u1[1] = u1[2] = u2[0] = u2[1] = u2[2] = 0.0;
  nod1 = 0L;
  nod2 = 0L;
  dcos = dsin = snorml = 0.0;    
}



axisrotrec::~axisrotrec()
{
}



/**
  The function reads one record for rotation of selected nodes about the given axis.
  The axis is given by two points A, B and rotation angle is given by directional cosine and sine. 
  The function is used in the growing mechanical problems
  in the case that the initial displacements of new active parts of structure should be calculated. Fo each new independent 
  part of structure that is being activated, one record consisting of time, axis definition, selected nodes, and directional 
  cosines may be specified. The selected nodes are expected to be at the surface of the new active part of structure and 
  the vector of their prescribed initial displacements will obtained by the difference between the original position and 
  position obtained by rotation of selected nodes about the given axis.

  @param in - pointer to the opened XFILE structure

  @return The function does not return anything but it stores the input data in the data members of the class.

  Created by Tomas Koudelka, 07.2016
*/
void axisrotrec::read(XFILE *in)
{
  long i;
  vector tmp(ASTCKVEC(3));

  read_prep(in); // reading of rotation parameters
  seln.read(in); // reading of selected nodes
  switch (dax)
  {
    case ptx2:
      break;
    case nod_pt:
      Gtm->give_node_coord(nod1, tmp);
      copyv(tmp, a);
      break;
    case nodx2:
      Gtm->give_node_coord(nod1, tmp);
      copyv(tmp, a);
      Gtm->give_node_coord(nod2, tmp);
      copyv(tmp, b);
      break;
    default:
      print_err("unknown type of axis definition is required", __FILE__, __LINE__, __func__);
  }

  for(i=0; i<3; i++)
    l[i] = b[i]-a[i];
  scprd(l, l, snorml);
  if (snorml < Mp->zero)
    print_err("definition points of axis are too close to each other\n,"
              " square of directional vector norm %le < %le", __FILE__, __LINE__, __func__, snorml, Mp->zero);
  norml = sqrt(snorml);
}



/**
  The function reads parameters of rotation from the opened text file.

  @param in - pointer to the opened XFILE structure

  @return The function does not return anything but it stores the input data in the data members of the class.

  Created by Tomas Koudelka, 07.2016
*/
void axisrotrec::read_prep(XFILE *in)
{
  xfscanf(in, "%k%le%k%le", "start_time", &stime, "end_time", &etime);
  xfscanf(in, "%k%m", "def_axis", &axisdeftype_kwdset, &dax);
  switch (dax)
  {
    case ptx2:
      xfscanf(in, "%k%le%le%le", "coord_a", a, a+1, a+2);  
      xfscanf(in, "%k%le%le%le", "coord_b", b, b+1, b+2);
      break;
    case nod_pt:
      xfscanf(in, "%k%ld", "node_id1", &nod1);
      nod1--;
      xfscanf(in, "%k%le%le%le", "coord_b", b, b+1, b+2);
      break;
    case nodx2:
      xfscanf(in, "%k%ld", "node_id1", &nod1);
      nod1--;
      xfscanf(in, "%k%ld", "node_id2", &nod2);
      nod2--;
      break;
    default:
      print_err("unknown type of axis definition is required", __FILE__, __LINE__, __func__);
  }
  xfscanf(in, "%k%le%k%le", "dcos", &dcos, "dsin", &dsin);
  xfscanf(in, "%k%m", "apply_at_x", &answertype_kwdset, dap);
  xfscanf(in, "%k%m", "apply_at_y", &answertype_kwdset, dap+1);
  if (Gtm->give_maxdimension() > 2)
    xfscanf(in, "%k%m", "apply_at_z", &answertype_kwdset, dap+2);
}



/**
  The function prints a record of input data to the opened text file.

  @param out - pointer to the opened text ouput file
 
  @return The function does not return anything.

  Created by Tomas Koudelka, 07.2016
*/
void axisrotrec::print(FILE *out)
{
  fprintf(out, "%le %le\n", stime, etime);
  fprintf(out, "%d\n", dax);
  switch(dax)
  {
    case ptx2:
      fprintf(out, "%le %le %le\n", a[0], a[1], a[2]);
      fprintf(out, "%le %le %le\n", b[0], b[1], b[2]);
      break;
    case nod_pt:     
      fprintf(out, "%ld %le %le %le\n", nod1+1, b[0], b[1], b[2]);
      break;
    case nodx2:     
      fprintf(out, "%ld %ld\n", nod1+1, nod2+1);
      break;
    default:
      print_err("unknown type of axis definition is required", __FILE__, __LINE__, __func__);
  }
  fprintf(out, "%le %le\n", dcos, dsin);
  fprintf(out, "%d ", dap[0]);
  fprintf(out, "%d ", dap[1]);
  if (Gtm->give_maxdimension() > 2)
    fprintf(out, "%d\n", dap[2]);
  else
    fprintf(out, "\n");

  seln.print(out);
}


/**
  The function stores actual displacement %vector at nodes used for the definition of 
  active rotation axis accoridng to the actual time. If the axis is being defined with 
  the help of 2 points then nothing is performed.

  @param lcid - load case id
  @param time - actual time

  @return The function does not return anything but it stores retrievd displacement values
          to the class arrays u1 and u2.

  Created by Tomas Koudelka 08.2016
*/
void axisrotrec::store_displ_def_node(long lcid, double time)
{
  long aux, ndofn;
  vector d;

  // test if the rotation axis is active
  if ((time < stime) || (time > etime))
    return;

  switch (dax)
  {
    case ptx2:
      break;
    case nod_pt:
      // store displacement vector from nod1
      memset(u1, 0, sizeof(*u1)*3);
      ndofn = Mt->give_ndofn(nod1);
      if (ndofn < 3)
        reallocv(RSTCKVEC(3, d));
      else
        reallocv(RSTCKVEC(ndofn, d));
      noddispl(lcid, d.a, nod1);
      d.n = 3;
      copyv(d, u1);
      break;
    case nodx2:
      // store displacement vectors from nod1 and nod2
      memset(u1, 0, sizeof(*u1)*3);
      memset(u2, 0, sizeof(*u2)*3);
      ndofn = Mt->give_ndofn(nod1);
      if (ndofn < 3)
        reallocv(RSTCKVEC(3, d));
      else
        reallocv(RSTCKVEC(ndofn, d));
      noddispl(lcid, d.a, nod1);
      aux = d.n;
      d.n = 3;
      copyv(d, u1);
      d.n = aux;
      noddispl(lcid, d.a, nod2);
      d.n = 3;
      copyv(d, u2);
      break;
    default:
      print_err("unknown type of axis definition is required", __FILE__, __LINE__, __func__);
  }
  return;
}



/**
  The function computes prescribed initial displacement of one node whose spatial coordinates 
  are given in %vector u by rotation of this node about the given axis where the rotation is given 
  by directional cosine/sine. Resulting displacements are stored at %vector ur.
  The function is intended for the usage at growing mechanical problems.

  @param u - coordinates of one selected node U stored in %vector (input)
  @param ur - resulting prescribed initial displacement %vector (ouput)

  @return The function returns resulting displacement vector at argumnet ur.

  Created by Tomas Koudelka, 07.2016
*/
void axisrotrec::compute_rotcoord(vector &u, vector &ur)
{
  vector r(ASTCKVEC(3));
  vector uid(ASTCKVEC(3));
  vector d(ASTCKVEC(3));
  matrix t(ASTCKMAT(3,3));
  double ldotc, normcu;
  

  // origin C of local coordinate system (x',y',z') where 
  // C := (AU)->(AB)
  vector c(ASTCKVEC(3));
  
  // x' <=> vector of the given axis lx=l=(AB), 
  // y' <=> vector (CU) = ly
  // z' is the cross product of l and (CU): lz = l x (CU)
  // base vectors of local coordinate system axis-node
  vector ly(ASTCKVEC(3));
  vector lz(ASTCKVEC(3));

  // vector AU  
  r[0] = u[0]-a[0];
  r[1] = u[1]-a[1];
  r[2] = u[2]-a[2];
 
  // Calculation of the origin C of the local coordinate system is obtained by projection of AU onto AB
  scprd(l, r, ldotc);
  c[0] = ldotc/snorml*l[0] + a[0];
  c[1] = ldotc/snorml*l[1] + a[1];
  c[2] = ldotc/snorml*l[2] + a[2];

  // Assemble transformation matrix
  // base vector of local axis x'
  t(0,0) = l[0]/norml; 
  t(0,1) = l[1]/norml; 
  t(0,2) = l[2]/norml; 
  // base vector of local axis y'
  ly[0] = u[0]-c[0];
  ly[1] = u[1]-c[1];
  ly[2] = u[2]-c[2];
  normcu = normv(ly);
  cmulv(1.0/normcu, ly);
  t(1,0) = ly[0]; 
  t(1,1) = ly[1]; 
  t(1,2) = ly[2]; 
  // base vector of local axis z'
  crprd(l, ly, lz);
  cmulv(1.0/normv(lz), lz);
  t(2,0) = lz[0]; 
  t(2,1) = lz[1]; 
  t(2,2) = lz[2]; 

  // Components of prescribed initial displacement vector at the local coord. system
  // u_{id} = u' - u'_r where u'={0, |CU|, 0}^T and u'_r={0, dcos.|CU|, dsin.|CU|}^T
  uid[0] = 0.0;
  uid[1] = (dcos-1.0)*normcu;
  uid[2] = dsin*normcu;
 
  // Transformation of prescribed initial displacement vector to the global coordinate system
  mtxv(t, uid, ur);

  switch (dax)
  {
    case ptx2:
      break;
    case nod_pt:
      // add displacement vector from nod1 to the initial displacements due to rotation
      addv(ur.a, u1, ur.a, 3);
      break;
    case nodx2:
      // add average displacement vector from nod1 and nod2 to the initial displacements due to rotation
      addv(u1, u2, d.a, 3);
      cmulv(0.5, d);
      addv(ur, d, ur);
      break;
    default:
      print_err("unknown type of axis definition is required", __FILE__, __LINE__, __func__);
  }
}
