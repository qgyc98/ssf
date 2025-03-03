#ifndef ELEMTOOLS_H
#define ELEMTOOLS_H

#include "vector.h"

/// function generates outer normal %vector to required element end
void bar_normal_vectors (long enid,vector &n);

///  function constructs the outer unit normal %vector to an element edge
void triangle_normal_vectors (long edgeid,vector &x,vector &y,vector &n);

///  function constructs the outer unit normal %vector to an element edge
void quadlin_normal_vectors (long edgeid,vector &x,vector &y,vector &n);

///  function computes lengths of the quadrilateral edges
double quadlin_edge_length (long edgeid,vector &x,vector &y);

///  function computes areas of the axisymmetric quadrilateral surfaces
double quadlinaxisym_surface_area (long surfid,vector &x,vector &y);

///  function computes volume of an element
double tetrahedra_volume (long eid,vector &x,vector &y,vector &z);

///  function constructs the outer unit normal %vector to an element surface
void tetrahedra_normal_vectors (long surfid,vector &x,vector &y,vector &z,vector &n);

///  function computes areas of the tetrahedron surfaces
double tetrahedra_surface_areas (long surfid,vector &x,vector &y,vector &z);

///  function computes volume of an element
double hexahedra_volume (long eid,vector &x,vector &y,vector &z);

///  function constructs the outer unit normal %vector to an element surface
void hexahedra_normal_vectors (long surfid,vector &x,vector &y,vector &z,vector &n);

///  function computes areas of the hexahedron surfaces
double hexahedra_surface_areas (long surfid,vector &x,vector &y,vector &z);
#endif
