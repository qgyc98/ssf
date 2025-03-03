#ifndef TENSORCNAME
#define TENSORCNAME

#include "alias.h"

/**
  Class that defines stress/strain tensor component names and
  index names of a general tensors for particular stress/strain states.

  Created by Tomas Koudelka, 10.2023
*/
class tensor_cname
{
 public:
  //
  // the second order tensor indices for particular stress/strain states
  //
  /// tensor indices for uniaxial stress/strain state
  static const char*tens_indices_bar[];
  /// tensor indices for stress/strain arrays on 2D beam
  static const char*tens_indices_plbeam[];
  /// tensor indices for stress/strain arrays on 3D beam
  static const char*tens_indices_spbeam[];
  /// tensor indices for plane stress/strain state
  static const char*tens_indices_plstr[];
  /// tensor indices for plane contact stress/strain state
  static const char*tens_indices_plcont[];
  /// tensor indices for stress/strain state in Kirchoff's plates
  static const char*tens_indices_platek[];
  /// tensor indices for stress/strain state in Mindlin's plates
  static const char*tens_indices_plates[];
  /// tensor indices for axisymmetric stress/strain state
  static const char*tens_indices_axisym[];
  /// tensor indices for stress/strain state in shells
  static const char*tens_indices_shell[];
  /// tensor indices for fully 3D stress/strain state
  static const char*tens_indices_space[];

  //
  // strain tensor component names for particular strain states, Voigt notation is considered
  //
  /// strain tensor component names for uniaxial strain state
  static const char*strain_comp_names_bar[];
  /// strain tensor component names of strains for 2D beams
  static const char*strain_comp_names_plbeam[];
  /// strain tensor component names of strains for 3D beams
  static const char*strain_comp_names_spbeam[];
  /// strain tensor component names for plane strain state
  static const char*strain_comp_names_plstr[];
  /// strain tensor component names for plane constact strain state
  static const char*strain_comp_names_plcont[];
  /// curvature component names for Krichof's plate
  static const char*strain_comp_names_platek[];
  /// curvature nad shear strain component names for Mindlin's plate
  static const char*strain_comp_names_plates[];
  /// strain tensor component names for axisymmetric strain state
  static const char*strain_comp_names_axisym[];
  /// strain tensor component names for strain state on shells
  static const char*strain_comp_names_shell[];
  /// strain tensor component names for fully 3D strain state
  static const char*strain_comp_names_space[];

  //
  // stress tensor component names for particular stress states, Voigt notation is considered
  //
  /// stress tensor component names for uniaxial stress state
  static const char*stress_comp_names_bar[];
  /// strain tensor component names of strains for 2D beams
  static const char*stress_comp_names_plbeam[];
  /// strain tensor component names of strains for 3D beams
  static const char*stress_comp_names_spbeam[];
  /// stress tensor component names for plane stress state
  static const char*stress_comp_names_plstr[];
  /// stress tensor component names for plane constact stress state
  static const char*stress_comp_names_plcont[];
  /// stress resultant component names for Krichof's plate
  static const char*stress_comp_names_platek[];
  /// stress resultant component names for Mindlin's plate
  static const char*stress_comp_names_plates[];
  /// stress tensor component names for axisymmetric stress state
  static const char*stress_comp_names_axisym[];
  /// stress resultant component names on shells
  static const char*stress_comp_names_shell[];
  /// stress tensor component names for fully 3D stress state
  static const char*stress_comp_names_space[];
  
  /// returns pointer to string with indices for the given second order tensor component
  static const char* tens_indstr(strastrestate ssst, long compid);
  /// returns pointer to string with label for the given second order stress tensor component
  static const char* stress_cmpstr(strastrestate ssst, long compid);
  /// returns pointer to string with label for the given second order strain tensor component
  static const char* strain_cmpstr(strastrestate ssst, long compid);
};

#endif
