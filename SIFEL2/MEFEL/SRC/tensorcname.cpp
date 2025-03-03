#include "tensorcname.h"
#include "iotools.h"

#include <stdlib.h>

// strain component names for particular strain states
const char*tensor_cname::strain_comp_names_bar[]    = {"eps_x"};
const char*tensor_cname::strain_comp_names_plbeam[] = {"u_l", "w_l", "phi"};
const char*tensor_cname::strain_comp_names_spbeam[] = {"u_l", "v_l", "w_l", "phi_xl", "phi_yl", "phi_zl"};
const char*tensor_cname::strain_comp_names_plstr[]  = {"eps_x", "eps_y", "gamma_xy", "eps_z"};
const char*tensor_cname::strain_comp_names_plcont[] = {"eps_n", "eps_tx", "eps_ty"};
const char*tensor_cname::strain_comp_names_platek[] = {"kappa_x", "kappa_y", "kappa_xy"};
const char*tensor_cname::strain_comp_names_plates[] = {"kappa_x", "kappa_y", "kappa_xy", "gamma_xz", "gamma_yz"};
const char*tensor_cname::strain_comp_names_axisym[] = {"eps_r", "eps_y", "eps_phi", "gamma_ry"};
const char*tensor_cname::strain_comp_names_shell[]  = {"eps_x", "eps_y", "gamma_xy", "kappa_x", "kappa_y", "kappa_xy"};
const char*tensor_cname::strain_comp_names_space[]  = {"eps_x", "eps_y", "eps_z", "gamma_yz", "gamma_xz", "gamma_xy"};

// stress component names for particular stress states
const char*tensor_cname::stress_comp_names_bar[]    = {"sig_x"};
const char*tensor_cname::stress_comp_names_plbeam[] = {"N", "V", "M"};
const char*tensor_cname::stress_comp_names_spbeam[] = {"N", "V_yl", "V_zl", "M_xl", "M_yl", "M_zl"};
const char*tensor_cname::stress_comp_names_plstr[]  = {"sig_x", "sig_y", "tau_xy", "sig_z"};
const char*tensor_cname::stress_comp_names_plcont[] = {"sig_n", "tau_tx", "tau_ty"};
const char*tensor_cname::stress_comp_names_platek[] = {"m_x", "m_y", "m_xy"};
const char*tensor_cname::stress_comp_names_plates[] = {"m_x", "m_y", "m_xy", "q_xz", "q_yz"};
const char*tensor_cname::stress_comp_names_axisym[] = {"sig_r", "sig_y", "sig_phi", "tau_ry"};
const char*tensor_cname::stress_comp_names_shell[]  = {"n_x", "n_y", "n_xy", "m_x", "m_y", "m_xy"};
const char*tensor_cname::stress_comp_names_space[]  = {"sig_x", "sig_y", "sig_z", "tau_yz", "tau_xz", "tau_xy"};

// tensor component indices for particular strain/stress states
const char*tensor_cname::tens_indices_bar[]    = {"x"};
const char*tensor_cname::tens_indices_plbeam[] = {"x", "z", "phi"};
const char*tensor_cname::tens_indices_spbeam[] = {"xl", "yl", "zl", "xl", "yl", "zl"};
const char*tensor_cname::tens_indices_plstr[]  = {"x", "y", "xy", "z"};
const char*tensor_cname::tens_indices_plcont[] = {"n", "tx", "ty"};
const char*tensor_cname::tens_indices_platek[] = {"x", "y", "xy"};
const char*tensor_cname::tens_indices_plates[] = {"x", "y", "xy", "xz", "yz"};
const char*tensor_cname::tens_indices_axisym[] = {"r", "y", "phi", "ry"};
const char*tensor_cname::tens_indices_shell[]  = {"x", "y", "xy", "x", "y", "xy"};
const char*tensor_cname::tens_indices_space[]  = {"x", "y", "z", "yz", "xz", "xy"};



/** 
  The function returns pointer to string with indices for the given second order tensor component.

  @param[in] ssst - tensor stress/strain state
  @param[in] compid - required tensor component, tensor storage in Voigt notation is considered

  @return The function returns pointer to string representation of the given index of required
          tensor component.

  Created by Tomas Koudelka, 10.2023
*/
const char*tensor_cname::tens_indstr(strastrestate ssst, long compid)
{
  const char *ret = NULL;
  
  switch (ssst){
    case bar:
      ret = tens_indices_bar[compid];
      break;
    case planestress:
    case planestrain:
      ret = tens_indices_plstr[compid];
      break;
    case planecontact:
      ret = tens_indices_plcont[compid];
      break;
    case platek:
      ret = tens_indices_platek[compid];
      break;
    case plates:
      ret = tens_indices_plates[compid];
      break;
    case axisymm:
      ret = tens_indices_axisym[compid];
      break;
    case shell:
      ret = tens_indices_shell[compid];
      break;
    case spacestress:
      ret = tens_indices_space[compid];
      break;
    default:
      print_err("unknown stress/strain state %d is required", __FILE__, __LINE__, __func__, int(ssst));
      abort();
  }
  return ret;
}



/** 
  The function returns pointer to string with label for the given strain tensor component in Voigt notation.

  @param[in] ssst - tensor strain state
  @param[in] compid - required tensor component, tensor storage in Voigt notation is considered

  @return The function returns pointer to string representation of the required strain tensor component in Voigt notation.

  Created by Tomas Koudelka, 10.2023
*/
const char*tensor_cname::strain_cmpstr(strastrestate ssst, long compid)
{
  const char *ret = NULL;
  
  switch (ssst){
    case bar:
      ret = strain_comp_names_bar[compid];
      break;
    case planestress:
    case planestrain:
      ret = strain_comp_names_plstr[compid];
      break;
    case planecontact:
      ret = strain_comp_names_plcont[compid];
      break;
    case platek:
      ret = strain_comp_names_platek[compid];
      break;
    case plates:
      ret = strain_comp_names_plates[compid];
      break;
    case axisymm:
      ret = strain_comp_names_axisym[compid];
      break;
    case shell:
      ret = strain_comp_names_shell[compid];
      break;
    case spacestress:
      ret = strain_comp_names_space[compid];
      break;
    default:
      print_err("unknown strain state %d is required", __FILE__, __LINE__, __func__, int(ssst));
      abort();
  }
  return ret;
}



/** 
  The function returns pointer to string with label for the given stress tensor component in Voigt notation.

  @param[in] ssst - tensor stress state
  @param[in] compid - required tensor component, tensor storage in Voigt notation is considered

  @return The function returns pointer to string representation of the required stress tensor component in Voigt notation.

  Created by Tomas Koudelka, 10.2023
*/
const char*tensor_cname::stress_cmpstr(strastrestate ssst, long compid)
{
  const char *ret = NULL;
  
  switch (ssst){
    case bar:
      ret = stress_comp_names_bar[compid];
      break;
    case planestress:
    case planestrain:
      ret = stress_comp_names_plstr[compid];
      break;
    case planecontact:
      ret = stress_comp_names_plcont[compid];
      break;
    case platek:
      ret = stress_comp_names_platek[compid];
      break;
    case plates:
      ret = stress_comp_names_plates[compid];
      break;
    case axisymm:
      ret = stress_comp_names_axisym[compid];
      break;
    case shell:
      ret = stress_comp_names_shell[compid];
      break;
    case spacestress:
      ret = stress_comp_names_space[compid];
      break;
    default:
      print_err("unknown stress state %d is required", __FILE__, __LINE__, __func__, int(ssst));
      abort();
  }
  return ret;
}
