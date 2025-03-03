#ifndef INPUTT_H
#define INPUTT_H

#include <stdio.h>
#include "descript.h"
#include "xfile.h"
#include "siftop.h"
#include "prepalias.h"
#include "list.h"
#include "aliast.h"
#include "entitybocon.h"
#include "intdomt.h"

/// main function for input from the preprocessor file
long inputt(XFILE *in, descript *d);

/// input section with names of topology, material and cross section files
long input_filest(XFILE *in, descript &d);

/// input of description of load cases
long input_lct(XFILE *in);

/// input topolgy using siftop class
long input_siftopt(XFILE *in, descript *d);

/// input material from the material file
long input_materialst(char *fname, descript *d);

/// input cross-section from the cross-section file
long input_crst(char *fname, descript *d);

/// input material from the preprocessor file section
long input_materialst(XFILE *in, descript *d);

/// input cross-section from the preprocessor file section 
long input_crst(XFILE *in, descript *d);

/// input of nodal properties
long input_nodpropt(XFILE *in);

/// input of element properties
long input_elempropt(XFILE *in);

/// input of hanging nodes from the separate file
long input_hang_nodest(XFILE *in);

/// assignment of boundary conditions at nodes
long input_nod_bocont(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of coupled dofs at nodes
long input_nod_coupl_dofst(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of time functions cotrolling dofs at nodes
long input_nod_dof_tfunct(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of cross-setions at nodes
long input_nod_crsect(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of initial conditions at nodes
long input_nod_initcondt(XFILE *in, const enumstr nodsects[], long nsect);

/// reading of initial conditions from file
long read_inicd_file(char *cndfname, long ndof);

/// assignment of sources of quantities at nodes
long input_nod_sourcet(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of nodal advection velocities at nodes
long input_nod_advect_vel(XFILE *in, const enumstr nodsects[], long nsect);

/// assignment of element types
long input_elem_typet(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment of element material types
long input_elem_matt(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment of element cross-section types
long input_elem_crsect(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment of sources of quantities on elements
long input_elem_sourcet(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment vertex (nodal) load to elements
long input_elem_vertbc(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment edge load to elements
long input_elem_edgebc(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment surface load to elements
long input_elem_surfbc(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment of time functions to elements
long input_elem_eltimefunct(XFILE *in, const enumstr elemsects[], long nsect);

/// assignment of integration domain descriptors
long input_elem_fluxint(XFILE *in, const enumstr elemsects[], long nsect);

gentity operator --(gentity &a, int);

/// function assembles bnodval object from the entitybocon object and array of bndoval indeces for the given element
long assemble_bnodvalt(selement &el, snode *nodes, long prop, entitybocon &ebc, list &bnvl, gentity ent, long minid, long *bnid, long *entid, long nentid);

/// function assmebles intdomt object according to detected property id on the given element
void assemble_int_domain(selement &el, snode *nodes, long prop, gentity ent, long *entid, long nentid, long idomid, intdomt &idom);

/// function assembles array of indeces of climatcond object for the given element
void assemble_bclimcond(selement &el, snode *nodes, long prop, list &bclimc, gentity ent, long *bnid, long *entid, long nentid);

/// function returns pointer to a new loadelt structure with BC created for the given entity bc, element and entity property id
loadelt *bc2loadelt(long eid, selement &el, snode *nodes, long prop, bocontypet bc, long **bnid, gentity ent, long *entid, long nentid);

/// returns the number of boundary objects 
long get_nbo(selement &el, gentity ent);

long merge_intdom(list *idoml, intdomt *idom);
#endif
