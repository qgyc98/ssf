#ifndef GLOBPREPT_H
#define GLOBPREPT_H

#include "vector.h"
#include "siftop.h"
#include "list.h"
#include "aliast.h"
#include "globalt.h"
#include "hangnodet.h"
#include "bocont.h"
#include "sourcet.h"
#include "loadelt.h"
#include "dbcrst.h"
#include "dbmatt.h"
#include "entitybocon.h"


#ifndef EXTERN
 #define EXTERN extern
#endif

EXTERN char        *Inicdf;  ///< pointer to the initial condition file name (must not be deleted)
EXTERN FILE        *Log;     ///< Log file of mechprep
EXTERN siftop      *Top;     ///< Sifel topology
EXTERN dbmatt      *Dbmatt;  ///< database of materials and their parameters
EXTERN dbcrst      *Dbcrst;  ///< database of cross-sections and their parameters

EXTERN gfunct      *Tft;     ///< array of time functions used in time dependent mechanical problems for each load case

EXTERN long         Numhn;   ///< number of hanging nodes

// Following arrays have size equal to number of nodes (first index has to be node number)
EXTERN hangnodet  **Nod_hang;      ///< connectivity of hanging nodes
EXTERN long       **Nod_bocon;     ///< boundary conditions at nodes (prescribed values of unknowns) Nod_bc[i][j] returns id of pvalt for i-th node and j-th lcid
EXTERN list         Nod_bclst;     ///< list of particular prescribed values of uknowns. It is referenced by Nod_bocon
EXTERN long       **Nod_ccn;       ///< prescribed coupled dofs at nodes
EXTERN crsectypet  *Nod_cst;       ///< cross-section types at nodes
EXTERN long        *Nod_csti;      ///< cross-section indeces at nodes
EXTERN long        *Nod_cstdbi;    ///< cross-section type indeces of nodes in cross-section database
EXTERN double     **Nod_inicd;     ///< initial conditions at nodes, Nod_inicd[nod_id][lcid] = value of the initial condition for the given quantity (lcid) and given node (nod_id)
EXTERN long       **Nod_sourcet;   ///< nodal sources of transported quantities, Nod_sourcet[nod_id][lcid] = index of quantitiy source in the Nod_srclst
EXTERN list         Nod_srclst;    ///< list of particular nodal sources of quantities, indeces from Nod_sourcet refer to this list
EXTERN long        *Nod_advect_nc; ///< The number of advection velocity components at nodes, Nod_advect_nc[nod_id] = number of velocity components
EXTERN double     **Nod_advect;    ///< Nodal vectors of advection velocities, Nod_advect[nod_id][j] = j-th advection velocity component at node nod_id

// Following arrays have size equal to number of elements (first index has to be element number)
EXTERN elemtypet     *El_type;   ///< type of elements
EXTERN long          *El_nmat;   ///< number of material types assigned to particular elements
EXTERN mattypet     **El_mattype;///< material type of elements
EXTERN long         **El_matid;  ///< material indeces of elements
EXTERN long         **El_matdbi; ///< material type indeces of elements in material database
EXTERN crsectypet    *El_cst;    ///< cross-section of elements
EXTERN long          *El_csti;   ///< cross-section indeces of elements
EXTERN long          *El_cstdbi; ///< cross-section type indeces of elements in cross-section database
EXTERN long         **El_sourcet;///< sources of transported quantities on elements, El_sourcet[elem_id][lcid] = index of quantitiy source in the El_srclst
EXTERN list           El_srclst; ///< list of particular sources of quantities on elements, indeces from El_sourcet refer to this list

EXTERN loadelt     ***El_loadt;     ///< element boundary condition, El_loadt[lcid][el_id] = pointer to loadelt structure
EXTERN long         **El_loadtln;   ///< lines where the vertex|edge|surface El_loadt has been assigned firstly, El_loadtln[lcid][el_id] first line with assignment of vertex|edge|surface boundary condition
EXTERN long         **El_loadtcol;  ///< columns where the vertex|edge|surface El_loadt has been assigned firstly, El_loadtcol[lcid][el_id] first column with assignment of vertex|edge|surface boundary condition
EXTERN list          *El_nv_lst;    ///< El_nv_lst[lcid]    = list of bnodvalt objects for nodal values whose indeces are used in El_loadt for the load case lcid
EXTERN list          *El_trc_lst;   ///< El_trc_lst[lcid]   = list of bnodvalt objects for transmission coefficients whose indeces are used in El_loadt for the load case lcid 
EXTERN list          *El_trr_lst;   ///< El_trr_lst[lcid]   = list of bnodvalt objects for radiation coefficients whose indeces are used in El_loadt for the load case lcid  
EXTERN list          *El_cc_lst;    ///< El_cc_lst[lcid]    = list of climatcond objects whose indeces are used in El_loadt for the load case lcid 
EXTERN list          *El_gcc_lst;   ///< El_gcc_lst[lcid]   = list of climatcond2 objects whose indeces are used in El_loadt for the load case lcid 
EXTERN list          *El_ccf_lst;   ///< El_ccf_lst[lcid]   = list of files with climatcond objects whose indeces are used in El_loadt for the load case lcid 
EXTERN list          *El_gccf_lst;  ///< El_gccf_lst[lcid]  = list of files with climatcond2 objects whose indeces are used in El_loadt for the load case lcid 
EXTERN list          *El_trcc_lst;  ///< El_trcc_lst[lcid]  = list of types of climatcond objects whose indeces are used in El_loadt for the load case lcid 
EXTERN list          *El_gtrcc_lst; ///< El_gtrcc_lst[lcid] = list of types of climatcond2 objects whose indeces are used in El_loadt for the load case lcid 
EXTERN list        ***El_fluxint;   ///< El_fluxint[lcid][el_id] = pointer to the list of intdom structure with description of domain integral indeces
EXTERN long           Nidomid;      ///< the total number of flux integration domains (maximum integration domain id)
EXTERN long          *El_tfunc;     ///< indeces of time fucntions for switching elements on/off

// Auxiliary variables
EXTERN long   Check_unused_prop; ///< Flag for checking of unused properties

#endif
