#ifndef OUTPUT_H
#define OUTPUT_H
#include <stdio.h>
#include "descrip.h"

/// main output function
long output(FILE *out, descrip *d);

/// output of nodes
long wr_nodes(FILE *out);

/// output of nodes for layered problems
long wr_lay_nodes(FILE *out);

/// output of boundary conditions
long wr_bocon(FILE *out);

/// output of elements
long wr_elements(FILE *out, descrip *d);

/// output of elements for layered problems
long wr_lay_elements(FILE *out, descrip *d);

/// output of global node numbers for paralell computing
long wr_globnodnum(FILE *out, descrip *d);

/// output of periodic boundary condition data
long wr_periodbc(FILE *out);

/// output of integration point section (already not used)
long wr_intpoints(FILE *out, long nprop);

/// output of materials
long wr_materials(FILE *out, descrip *d);

/// output auxiliary point section
long wr_auxpoint(FILE *out);

/// output of section with local coordinate systems on the elements
long wr_ellcsys(FILE *out);

/// output of cross-sections
long wr_crsecs(FILE *out, descrip *d);

/// load controling function
long wr_load(FILE *out);

/// output of nodal load
long wr_loadn(FILE *out, long nlc);

/// output of element load
long wr_loadel(FILE *out, long nlc);

/// output of nodal dynamical load
long wr_dloadn(FILE *out, long nlc);

/// output of element dynamical load
long wr_dloadel(FILE *out, long nlc);

/// output of prescribed displacements
long wr_prescdisp(FILE *out, long nlc);

/// output of macro strain/stress components for homogenization problems
long wr_macrostrastre(FILE *out, long nlc);

/// output of dynamical prescribed displacements
long wr_dprescdisp(FILE *out, long nlc);

/// output of temperature load
long wr_tempload(FILE *out, long nlc);

/// output of initial conditions
long wr_initcond(FILE *out, long nlc);

/// output of prescribed initial displacements by rotation
long wr_rotinidispl(FILE *out, long nlc);

/// output of element eigenstrains
long wr_eigenstrains(FILE *out);

/// returns the total number of spring supports at nodes
long get_nspring();

/// tranformation of periodic bc to the hanging nodes for output 
void transf_periodbc_hangnod();
#endif
