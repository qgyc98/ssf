#ifndef GLOBPREP_H
#define GLOBPREP_H

#include "vector.h"
#include "siftop.h"
#include "list.h"
#include "alias.h"
#include "global.h"
#include "loadel.h"
#include "dloadel.h"
#include "hangnode.h"
#include "bocon.h"
#include "inicd.h"
#include "dbcrs.h"
#include "dbmat.h"
#include "pointset.h"
#include "tempload.h"
#include "axisrotrec.h"
#include "msmap.h"

#ifndef EXTERN
 #define EXTERN extern
#endif

EXTERN FILE      *Log;     ///< Log file of mechprep
EXTERN siftop    *Top;     ///< Sifel topology
EXTERN dbmat     *Dbmat;   ///< database of materials and their parameters
EXTERN dbcrs     *Dbcrs;   ///< database of cross-sections and their parameters

EXTERN long       Nlc;     ///< number of loading cases
EXTERN long      *Nslc;    ///< number of loading cases
EXTERN long      *Nslc_cum;///< cumulative number of subload cases, Nslc_cum[i] = sum(0,j<i)(Nslc[j])
EXTERN long       Tnslc;   ///< total number of all subload cases
EXTERN long       Nlay;    ///< number of layers for layered problems
EXTERN long       Numpd;   ///< total number of values of prescribed displacements (both static and time dependent)
EXTERN long       Numhn;   ///< number of hanging nodes
EXTERN long       Ntempl;  ///< indicator of temprature load
EXTERN gfunct    *Tf;      ///< array of time functions used in time dependent mechanical problems for each load case
EXTERN dynload    Tdload;  ///< type of time dependent load
EXTERN long      *Tlt;     /**< array of temperature load type for every load case
                                usage: Tlt[i] = type of temperature load case for subloadces with cumulative number i */
EXTERN long      *Npd;     /**< array with numbers of prescribed displacements used in growing mechanical problems
                                usage: Npd[i] = number of prescribed displacements in subloadces with cumulative number i */
EXTERN double   **Spd;     /**< array with values of prescribed displacements used for particular subloadcases in growing mechanical problems
                                usage: Spd[i][j]= value of j-th prescribed displacement for subloadcase with cumulative number i */
EXTERN long       Nmstrc;  ///< number of macro stress components in homogenization problem
EXTERN double   **Mstrc;   ///< array of macro stress components for particular load cases Mstrc[lcid][comp_id]
EXTERN strastre **Mstrct;  ///< array of macro component types (strain/stress) for particular load cases Mstrct[lcid][comp_id], it is intended for Mp->homog == 9
EXTERN long       Nperbc;  ///< total number of nodes with the prescribed periodic boundary condition.

// Following arrays have size equal to number of nodes (first index has to be node number)
EXTERN hangnode   **Nod_hang;    ///< connectivity of hanging nodes
EXTERN long        *Nod_ndof;    ///< number of degrees of freedom (dofs) at nodes
EXTERN bocon      **Nod_bocon;   ///< boundary conditions at nodes (prescribed values of unknowns)
EXTERN long       **Nod_ccn;     ///< prescribed coupled dofs at nodes
EXTERN crsectype   *Nod_cst;     ///< cross-section types at nodes
EXTERN long        *Nod_csti;    ///< cross-section indeces at nodes
EXTERN long        *Nod_cstdbi;  ///< cross-section type indeces of nodes in cross-section database
EXTERN vector     **Nod_lcs;     ///< base vectors of local coordinate systems at nodes
EXTERN double    ***Nod_load;    ///< load components at nodes, Nod_load[nod_id][lcid][comp_id]
EXTERN gfunct    ***Nod_tdload;  ///< time dependent load components at nodes, Nod_tdload[nod_id][lcid][comp_id] = gfunct for given load component
EXTERN inicd     ***Nod_inicd;   ///< initial conditions at nodes, Nod_inicd[nod_id][lcid] = pointer to initial condition structure
EXTERN tempload  ***Nod_temper;  ///< temperature load at nodes, Nod_temper[nod_id][lcid] = pointer to the temperature load
EXTERN long       **Nod_nsprmat;    ///< Nod_nsprmat[nod_id][dir]         = number of materials for spring support at node nod_id and direction dir
EXTERN mattype   ***Nod_sprmattype; ///< Nod_sprmattype[nod_id][dir][idm] = idm-th material type for spring support at node nod_id and direction dir
EXTERN long      ***Nod_sprmatid;   ///< Nod_sprmatid[nod_id][dir][idm]   = idm-th material index for spring support at node nod_id and direction dir
EXTERN long      ***Nod_sprmatdbi;  ///< Nod_sprmatdbi[nod_id][dir][idm]  = index of idm-th material type in dbmat for spring support at node nod_id and direction dir
EXTERN long         Nod_rot_ipd_num;///< The number of axis specified for the prescribed intial displacements (it is the length of array Nod_rot_ipd_axis)
EXTERN axisrotrec  *Nod_rot_ipd_axis; ///< array of records for the assignment of initial prescribed displacements by rotation about axis
EXTERN msmap       *Nod_periodbc;   ///< Nod_periodbbc[nod_id] = structure with master/slave node data for periodic boundary condition

// Following arrays have size equal to number of elements (first index has to be element number)
EXTERN elemtype      *El_type;     ///< type of elements
EXTERN strastrestate *El_ssst;     ///< strain/stress state of elements
EXTERN long          *El_nmat;     ///< number of material types assigned to particular elements
EXTERN mattype      **El_mattype;  ///< material type of elements
EXTERN long         **El_matid;    ///< material indeces of elements
EXTERN long         **El_matdbi;   ///< material type indeces of elements in material database
EXTERN crsectype     *El_cst;      ///< cross-section of elements
EXTERN long          *El_csti;     ///< cross-section indeces of elements
EXTERN long          *El_cstdbi;   ///< cross-section type indeces of elements in cross-section database
EXTERN vector       **El_lcs;      ///< base vectors of local coordinate systems at elements
EXTERN loadel      ***El_load;     ///< element load, El_load[el_id][lcid] = pointer to loadel structure
EXTERN long        ***El_loadln;   ///< lines where the edge|surface|volume El_load has been assigned firstly, El_loadln[el_id][lcid][0|1|2] first line with assignment of edge|surface|volume load
EXTERN long        ***El_loadcol;  ///< columns where the edge|surface|volume El_load has been assigned firstly, El_loadln[el_id][lcid][0|1|2] first column with assignment of edge|surface|volume load
EXTERN dloadel     ***El_tdload;   ///< element time dependent load, El_tdload[el_id][lcid] = pointer to dloadel structure
EXTERN long        ***El_tdloadln; ///< lines where the time dependent edge|surface|volume El_tdload has been assigned firstly, El_tdloadln[el_id][lcid][0|1|2] first line with assignment of edge|surface|volume load
EXTERN long        ***El_tdloadcol;///< columns where the time dependent edge|surface|volume El_tdload has been assigned firstly, El_loadln[el_id][lcid][0|1|2] first column with assignment of edge|surface|volume load
EXTERN long          *El_tfunc;    ///< indeces of time fucntions for switching elements on/off
EXTERN long         **El_eigstr;   ///< indeces of gfuncts for particular eigenstrain components on elements, El_eigstr[el_id][k] = index of gfunct in list El_eigstrgf_lst which defines the k-th eigenstrain component
EXTERN strastre      *El_eigstrt;  ///< type of assigned eigenquantity, i.e. eigenstrain or eigenstress
EXTERN list           El_eigstrgf_lst;///< list of gfuncts for eigenstrains on elements

// Auxiliary variables
EXTERN long   Check_unused_prop; ///< Flag for checking of unused properties

#endif
