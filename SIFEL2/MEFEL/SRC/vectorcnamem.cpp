#include "vectorcnamem.h"
#include "iotools.h"

#include <stdlib.h>


// names of a general %vector component indices for particular states
const char*vector_cnamem::vect_indices_bar[]    = {"x", "y", "z"};
const char*vector_cnamem::vect_indices_plbeam[] = {"x", "z", "ry"};
const char*vector_cnamem::vect_indices_spbeam[] = {"x", "y", "z", "rx", "ry", "ry"};
const char*vector_cnamem::vect_indices_plane[]  = {"x", "y"};
const char*vector_cnamem::vect_indices_plate[]  = {"z", "rx", "ry"};;
const char*vector_cnamem::vect_indices_shell[]  = {"x", "y", "z", "rx", "ry", "ry"};
const char*vector_cnamem::vect_indices_space[]  = {"x", "y", "z"};

// names of a displacement %vector indices for particular states
const char*vector_cnamem::displvec_comp_names_bar[]    = {"u", "v", "w"};
const char*vector_cnamem::displvec_comp_names_plbeam[] = {"u", "w", "phi_y"};
const char*vector_cnamem::displvec_comp_names_spbeam[] = {"u", "v", "w", "phi_x", "phi_y", "phi_z"};
const char*vector_cnamem::displvec_comp_names_plane[]  = {"u", "v"};
const char*vector_cnamem::displvec_comp_names_plate[]  = {"w", "phi_x", "phi_y"};
const char*vector_cnamem::displvec_comp_names_shell[]  = {"u", "v", "w", "phi_x", "phi_y", "phi_z"};
const char*vector_cnamem::displvec_comp_names_space[]  = {"u", "v", "w"};

// names of a load %vector indices for particular states
const char*vector_cnamem::reactvec_comp_names_bar[]    = {"R_x", "R_y", "R_z"};
const char*vector_cnamem::reactvec_comp_names_plbeam[] = {"R_x", "R_z", "Mr_y"};
const char*vector_cnamem::reactvec_comp_names_spbeam[] = {"R_x", "R_y", "R_z", "Mr_x", "Mr_y", "Mr_z"};
const char*vector_cnamem::reactvec_comp_names_plane[]  = {"R_x", "R_y"};
const char*vector_cnamem::reactvec_comp_names_plate[]  = {"R_z", "Mr_x", "Mr_y"};
const char*vector_cnamem::reactvec_comp_names_shell[]  = {"R_x", "R_y", "R_z", "Mr_x", "Mr_y", "Mr_z"};
const char*vector_cnamem::reactvec_comp_names_space[]  = {"R_x", "R_y", "R_z"};



/** 
  The function returns pointer to string with indices for the given vector component.

  @param[in] ssst - tensor stress/strain state
  @param[in] compid - required %vector component

  @return The function returns pointer to string representation of the required %vector component

  Created by Tomas Koudelka, 10.2023
*/
const char*vector_cnamem::vect_indstr(strastrestate ssst, long compid)
{
  const char *ret = NULL;
  
  switch (ssst){
    case bar:
      ret = vect_indices_bar[compid];
      break;
    case plbeam:
      ret = vect_indices_plbeam[compid];
      break;
    case spacebeam:
      ret = vect_indices_spbeam[compid];
      break;
    case planestress:
    case planestrain:
    case axisymm:
    case planecontact:
      ret = vect_indices_plane[compid];
      break;
    case platek:
    case plates:
      ret = vect_indices_plate[compid];
      break;
    case shell:
      ret = vect_indices_shell[compid];
      break;
    case spacestress:
      ret = vect_indices_space[compid];
      break;
    default:
      print_err("unknown stress/strain state %d is required", __FILE__, __LINE__, __func__, int(ssst));
      abort();
  }
  return ret;
}



/** 
  The function returns pointer to string with label for the given displacement vector component.

  @param[in] ssst - tensor stress state
  @param[in] compid - required displacement %vector component

  @return The function returns pointer to string representation of the required displacement %vector component.

  Created by Tomas Koudelka, 10.2023
*/
const char*vector_cnamem::displ_cmpstr(strastrestate ssst, long compid)
{
  const char *ret = NULL;
  
  switch (ssst){
    case bar:
      ret = displvec_comp_names_bar[compid];
      break;
    case plbeam:
      ret = displvec_comp_names_plbeam[compid];
      break;
    case spacebeam:
      ret = displvec_comp_names_spbeam[compid];
      break;
    case planestress:
    case planestrain:
    case planecontact:
    case axisymm:
      ret = displvec_comp_names_plane[compid];
      break;
    case plates:
    case platek:
      ret = displvec_comp_names_plate[compid];
      break;
    case shell:
      ret = displvec_comp_names_shell[compid];
      break;
    case spacestress:
      ret = displvec_comp_names_space[compid];
      break;
    default:
      print_err("unknown stress state %d is required", __FILE__, __LINE__, __func__, int(ssst));
      abort();
  }
  return ret;
}



/** 
  The function returns pointer to string with label for the given stress tensor component in Voigt notation.

  @param[in] ssst - tensor stress state
  @param[in] compid - required tensor component, tensor storage in Voigt notation is considered

  @return The function returns pointer to string representation of the required tensor component

  Created by Tomas Koudelka, 10.2023
*/
const char*vector_cnamem::react_cmpstr(strastrestate ssst, long compid)
{
  const char *ret = NULL;
  
  switch (ssst){
    case bar:
      ret = reactvec_comp_names_bar[compid];
      break;
    case plbeam:
      ret = reactvec_comp_names_plbeam[compid];
      break;
    case spacebeam:
      ret = reactvec_comp_names_spbeam[compid];
      break;
    case planestress:
    case planestrain:
    case planecontact:
    case axisymm:
      ret = reactvec_comp_names_plane[compid];
      break;
    case platek:
    case plates:
      ret = reactvec_comp_names_plate[compid];
      break;
    case shell:
      ret = reactvec_comp_names_shell[compid];
      break;
    case spacestress:
      ret = reactvec_comp_names_space[compid];
      break;
    default:
      print_err("unknown strain state %d is required", __FILE__, __LINE__, __func__, int(ssst));
      abort();
  }
  return ret;
}
