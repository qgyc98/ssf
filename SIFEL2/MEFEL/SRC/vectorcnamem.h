#ifndef VECTORCNAMEM
#define VECTORCNAMEM

#include "alias.h"


/**
  Class that defines nodal %vector component names for mechanical problems.
  It involves index names of a general %vectors and particular component names for 
  the displacement and reaction vectors that are defined with respect to the given 
  stress/strain state.

  Created by Tomas Koudelka, 10.2023
*/
class vector_cnamem
{
 public:
  //
  // General %vector indices for particular states
  //
  /// vector indices for uniaxial stress/strain state
  static const char*vect_indices_bar[];
  /// vector indices for arrays on 2D beam
  static const char*vect_indices_plbeam[];
  /// vector indices for arrays on 3D beam
  static const char*vect_indices_spbeam[];
  /// vector indices for plane stress/strain, plane contact and axisymmetric states
  static const char*vect_indices_plane[];
  /// vector indices for stress/strain state in Kirchoff's and Mindlin plates
  static const char*vect_indices_plate[];
  /// vector indices for shell state
  static const char*vect_indices_shell[];
  /// vector indices for fully 3D state
  static const char*vect_indices_space[];

  //
  // names of displacement %vector components for particular states
  //
  /// names of displacement %vector component for uniaxial stress/strain state
  static const char*displvec_comp_names_bar[];
  /// names of displacement %vector components for arrays on 2D beam
  static const char*displvec_comp_names_plbeam[];
  /// names of displacement %vector components for arrays on 3D beam
  static const char*displvec_comp_names_spbeam[];
  /// names of displacement %vector components for plane stress/strain, plane contact and axisymmetric states
  static const char*displvec_comp_names_plane[];
  /// names of displacement %vector components for stress/strain state in Kirchoff's and Mindlin plates
  static const char*displvec_comp_names_plate[];
  /// names of displacement %vector components for shell state
  static const char*displvec_comp_names_shell[];
  /// names of displacement %vector components for fully 3D state
  static const char*displvec_comp_names_space[];

  //
  // names of reaction %vector components for particular states
  //
  /// names of reaction %vector components for uniaxial stress/strain state
  static const char*reactvec_comp_names_bar[];
  /// names of reaction %vector components for arrays on 2D beam
  static const char*reactvec_comp_names_plbeam[];
  /// names of reaction %vector components for arrays on 3D beam
  static const char*reactvec_comp_names_spbeam[];
  /// names of reaction %vector components for plane stress/strain, plane contact and axisymmetric states
  static const char*reactvec_comp_names_plane[];
  /// names of reaction %vector components for stress/strain state in Kirchoff's and Mindlin plates
  static const char*reactvec_comp_names_plate[];
  /// names of reaction %vector components for shell state
  static const char*reactvec_comp_names_shell[];
  /// names of reaction %vector components for fully 3D state
  static const char*reactvec_comp_names_space[];

  /// returns pointer to string with indices for the given second order tensor component
  static const char*vect_indstr(strastrestate ssst, long compid);
  /// returns pointer to string with labels for the given displacement %vector component
  static const char*displ_cmpstr(strastrestate ssst, long compid);
  /// returns pointer to string with labels for the given load/reaction %vector component
  static const char*react_cmpstr(strastrestate ssst, long compid);
};

#endif
