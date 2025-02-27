#ifndef ORDERING_H
#define ORDERING_H

#include "vector.h"

void linbar_endpoints (long *nodes);
void quadbar_endpoints (long *nodes);

void lintriangle_edgnod (long *edgenod,long edg);
void quadtriangle_edgnod (long *edgenod,long edg);
void linquadrilat_edgnod (long *edgenod,long edg);
void quadquadrilat_edgnod (long *edgenod,long edg);
void cubicquadrilat_edgnod (long *edgenod,long edg);
void linhexahedral_edgnod (long *edgenod,long edg);
void quadhexahedral_edgnod (long *edgenod,long edg);
void lintetrahedral_edgnod (long *edgenod,long edg);
void quadtetrahedral_edgnod(long *edgenod,long edg);

void linhexahedral_surfnod (long *surfnod,long surf);
void quadhexahedral_surfnod (long *surfnod,long surf);
void lintetrahedral_surfnod (long *surfnod,long surf);
void quadtetrahedral_surfnod (long *surfnod,long surf);
void linwedge_surfnod (long *surfnod,long surf);

void nodcoord_bar (vector &xi);
void nodcoord_barq (vector &xi);
void nodcoord_planelt (vector &xi,vector &eta);
void nodcoord_planelq (vector &xi,vector &eta);
void nodcoord_planeqq (vector &xi,vector &eta);
void nodcoord_planecq (vector &xi,vector &eta);
void nodcoord_lintet (vector &xi,vector &eta,vector &zeta);
void nodcoord_quadtet (vector &xi,vector &eta,vector &zeta);
void nodcoord_linhex (vector &xi,vector &eta,vector &zeta);
void nodcoord_quadhex (vector &xi,vector &eta,vector &zeta);
void nodcoord_linwedge (vector &xi,vector &eta,vector &zeta);

void nodip_bar (long i,long n,ivector &ipnum);
void nodip_barq (long i,long n,ivector &ipnum);
void nodip_planelq (long i,long n,ivector &ipnum);
void nodip_planelt (long i,long n,ivector &ipnum);
void nodip_lintet (long i,ivector &ipnum);
void nodip_quadtet (long i,ivector &ipnum);
void nodip_linhex (long i,long n,ivector &ipnum);
void nodip_quadhex (long i,long n,ivector &ipnum);

#endif
