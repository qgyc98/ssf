#ifndef SEQUENT_H
#define SEQUENT_H

#include "alias.h"
#include <stdio.h>
#include <vector>
#include <set>

struct vector;
struct ivector;



/// creates list of adjacent integration points
void adjacip (void);

/// saves list of adjacent integration points
void save_adjacip();
/// restores list of adjacent integration points
long restore_adjacip();
/// saves list of adjacent integration points in a text file
void save_adjacip_txt(FILE *aux);
/// saves list of adjacent integration points in a binary file
void save_adjacip_bin(FILE *aux);
/// restores list of adjacent integration points in a text file
void restore_adjacip_txt(FILE *aux);
/// restores list of adjacent integration points in a binary file
void restore_adjacip_bin(FILE *aux);

/// returns estimated number of adjacent integration points
long give_max_adjacip();
/// reallocates memory for onedimensional array of long
void reallocate(long *&array, long req_size, long &capacity);
/// computes number of non-eliminated adjacent elements
void newnadjelel (long ipp,long *adjelelback,long *adjelel,long &nae,long &nnae);
/// creates list of adjacent elements
void newadjelel (long *adjelelback,long *adjelel,long &nae,long &nnae,long *treatedel,long &nte);
/// computes number of adjacent integration points closer than given limit
void in_dist (long ipp,vector &coord,double *ipc,long *adjelel,long nae,double rlim);
/// computes distances of adjacent integration points closer than given limit
void dist (long ipp,long max_nadjip,vector &coord,double *ipc,long *adjelel,long nae,double rlim);

/// creates list of adjacent integration points
void adjacip2 (void);
/// computes number of non-eliminated adjacent elements
void newnadjelel2 (const ivector &adjelel, ivector &adjelelback, long &nnae);
/// creates list of adjacent elements
void newadjelel2 (const ivector &adjelelback, std::set<long> &treatedel, ivector &adjelel, long &nnae);
/// computes number of adjacent integration points closer than given limit
void in_dist2 (long ipp, const vector &coord, const vector &ipc, double rlim, ivector &adjelel, std::vector<long> &ipl, std::vector<double> &ipd);
/// prints adjacent int. points to file
void print_adjacip(FILE *out);

#endif
