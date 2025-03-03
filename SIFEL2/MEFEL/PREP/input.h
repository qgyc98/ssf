#ifndef INPUT_H
#define INPUT_H

#include "descrip.h"
#include "xfile.h"
#include "prepalias.h"
#include "alias.h"
#include "hngen.h"
#include <stdio.h>

struct ivector;

/// main function for input from the preprocessor file
long input(XFILE *in, descrip *d);

/// reading of topology, material and cross section files
long input_files(XFILE *in, descrip &d);

/// input of description of load cases
long input_lc(XFILE *in);

/// input topolgy using siftop class
long input_siftop(XFILE *in, descrip *d);

/// input material from the material file
long input_materials(char *fname, descrip *d);

/// input material from the corresponding section of the preprocessor file
long input_materials(XFILE *in, descrip *d);

/// input cross-section from the cross-section file
long input_crs(char *fname, descrip *d);

/// input cross-section from the the corresponding section of the preprocessor file
long input_crs(XFILE *in, descrip *d);

/// input of nodal properties
long input_nodprop(XFILE *in);

/// input of element properties
long input_elemprop(XFILE *in);

/// input of hanging nodes from the separate file
long input_hang_nodes(XFILE *in);

/// assignment of dofs at nodes
long input_nod_ndof(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of boundary conditions at nodes
long input_nod_bocon(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of coupled dofs at nodes
long input_nod_coupl_dofs(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of time functions cotrolling dofs at nodes
long input_nod_dof_tfunc(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of cross-setions at nodes
long input_nod_crsec(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of spring supports at nodes
long input_nod_springs(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of local coordinate systems at nodes
long input_nod_lcs(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of nodal load
long input_nod_load(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of nodal time dependent load
long input_nod_tdload(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of initial conditions at nodes
long input_nod_initcond(XFILE *in, const enumstr nodsects[], long nsect);

/// assigment of initial prescribed displacements by rotation about the given axis
long input_nod_rotinidispl(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of temparture load at nodes
long input_nod_temper(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of periodic boundary conditions at nodes
long input_nod_periodic_bc(XFILE *in, const enumstr nodsects[], long nsect);

/// generation of hanging nodes 
long input_nod_genhn(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of element types
long input_elem_type(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment of element material types
long input_elem_mat(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment of element cross-section types
long input_elem_crsec(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment of local coordinate systems
long input_elem_lcs(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment of load defined at particular elements
long input_elem_load(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment edge load to elements
long input_elem_loadedge(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment surface load to elements
long input_elem_loadsurf(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment of volume load to elements
long input_elem_loadvol(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment of time dependent load defined at particular elements
long input_elem_tdload(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment time dependent edge load to elements
long input_elem_tdloadedge(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment time dependent surface load to elements
long input_elem_tdloadsurf(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment of time dependent volume load to elements
long input_elem_tdloadvol(XFILE *in, const enumstr elemsects[], long nsect);

/// asignment of eigenstrains to elements
long input_elem_eigstr(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment of time functions to elements
long input_elem_eltimefunc(XFILE *in, const enumstr elemsects[], long nsect);

gentity operator --(gentity &a, int);

gentity operator ++(gentity &a, int);

long give_lt_id(elloadtype tel);

void periodbc_log(long cmdid, XFILE *in, long aline, long acol, ivector &setsnodes, ivector &setmnodes, double rerr);

/// process, reduce and set array of master nodes on the given array of hanging nodes
void get_masternodes(long nn, hngen *hn);
#endif
