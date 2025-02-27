#ifndef OUTPUTT_H
#define OUTPUTT_H
#include <stdio.h>
#include "descript.h"

/// main output function
long outputt(FILE *out, descript *d);

/// output of nodes
long wr_nodest(FILE *out);

/// output of boundary conditions
long wr_bocont(FILE *out);

/// output of elements
long wr_elementst(FILE *out, descript *d);

/// output of global node numbers for paralell computing
long wr_globnodnumt(FILE *out);

/// output of materials
long wr_materialst(FILE *out, descript *d);

/// output of cross-sections
long wr_crsecst(FILE *out, descript *d);

/// load controling function
long wr_loadcaset(FILE *out);

/// output of prescribed quantities
long wr_prescquant(FILE *out);

/// output of sources at nodes and elements
long wr_sources(FILE *out, long nlc);

/// output of element load
long wr_loadelt(FILE *out, long nlc);

/// output of climatic conditions
long wr_climatcond(FILE *out, long nlc);

/// output of general climatic conditions
long wr_climatcond2(FILE *out, long nlc);

/// output of flux integration domain descriptors
long wr_fluxint(FILE *out);

/// output initial conditions
long wr_initcondt(FILE *out);

/// output of nodal advection velocity vectors
long wr_advect(FILE *out);
#endif
