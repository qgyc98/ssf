#ifndef AXISROTREC
#define AXISROTREC

#include <stdio.h>
#include "alias.h"
#include "selection.h"

struct XFILE;
struct vector;


/**
  The class is intended for the storage of one record for rotation of selected nodes about the given axis and calculation 
  of such a rotation. The axis is given by two points and rotation angle is given by directional cosine and sine. 
  This class is exploited in the growing mechanical problems in the case that the initial displacements of new active parts of 
  structure should be calculated. Fo each new independent part of structure that is being activated, one record consisting of time, 
  axis definition, selected nodes, and directional cosine/sine may be specified. The selected nodes are expected to be at the 
  surface of the new active part of structure and the vector of their prescribed initial displacements will obtained by the 
  difference between the original position and position obtained by rotation of selected nodes about the given axis.
*/
class axisrotrec
{
 public:
  /// start time at which the given rotation should applied
  double stime;
  /// end time at which the given rotation should applied
  double etime;
  /// type of axis definition (2 points, node+point, 2 nodes)
  axisdeftype dax;
  /// coordinates of the first axis point A
  double a[3];
  /// coordinates of the second axis point B
  double b[3];
  /// id of the first definition node of axis
  long nod1;
  /// displacement vector of the first definition node
  double u1[3];
  /// id of the second definition node of axis
  long nod2;
  /// displacement vector of the second definition node
  double u2[3];
  /// selection of nodes
  sel seln;
  /// directional cosine of the rotation angle
  double dcos;
  /// directional sine of the rotation angle
  double dsin;
  /// array of flags for directions where to apply the resulting prescribed initial displacements
  answertype dap[3];

  /// directional vector of axis l := B-A
  vector l;
  /// square of vector l norm
  double snorml;
  /// norm of %vector l
  double norml;

  axisrotrec();
  ~axisrotrec();
  /// reads input record (start_time, end_time, A, B, dcos, dsin and node selection) from file in and initializes l and snorml
  void read(XFILE *in);
  /// reads parameters of rotation (start_time, end_time, A, B, dcos, dsin) from the text file
  void read_prep(XFILE *in);
  /// prints input record to the file out
  void print(FILE *out);
  /// stores actual displacements at the active axis definition nodes
  void store_displ_def_node(long lcid, double time);
  /// computes new coordinates ur of point u rotated about the axis l
  void compute_rotcoord(vector &u, vector &ur);
};


#endif
