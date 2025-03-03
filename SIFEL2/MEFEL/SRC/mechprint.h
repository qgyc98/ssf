#ifndef MECHPRINT_H
#define MECHPRINT_H

#include "alias.h"
#include <stdio.h>
#include <set>


class gtopology;
class probdesc;
class mechtop;
class mechmat;
class siftop;


void print_init(long istep, const char *mode, long idn1=1, long ide1=1);
void print_step(long lcid, long istep, double lambda, double *fi);
void print_step(long lcid, long istep, double lambda, double *fi, double *fr);
void print_step_forced(long lcid, long istep, double lambda, double *fi);
void print_step_forced(long lcid, long istep, double lambda, double *fi, double *fr);
void print_flush();
void print_close();

void export_gid_mesh(FILE *out, long idn1, long ide1);
void export_gid_2dmesh(FILE *out, long icut, long idn1, long ide1);
void export_gid_gauss_pt(FILE *out, long ide1);
void write_gid_nodes(FILE *out, long idn1);
void write_gid_element(FILE *out, long i, long idn1, long ide1);
void write_gid_contact_element(FILE *out, long i, long idn1, long ide1);
void write_gid_2delement(FILE *out, long i, long id1, long nne, long icut, long di, long idn1, long ide1);
void write_gid_displ(FILE *out, long lcid, const char *desclcid);
void write_gid_nodscalar(FILE *out, strastre scal, long lcid, long dir, const char *desclcid);
void write_gid_nodvector(FILE *out, strastre q, long lcid, long sid, const char *desclcid);
void write_gid_nodtensor(FILE *out, strastre q, long lcid, long sid, const char *desclcid);
void write_gid_nforces(FILE *out, long lcid, const char *desclcid, const char *veclabel, double *ifor, bool print_react);
void write_gid_elemscalar(FILE *out, strastre scal, long lcid, long dir, const char *desclcid);
void write_gid_elem_type_scalar(FILE *out, strastre scal, long lcid, long dir, const char *desclcid, elemtype te);
void write_gid_elemvector(FILE *out, strastre q, long lcid, long sid, const char *desclcid);
void write_gid_elem_type_vector(FILE *out, strastre q, long lcid, long sid, const char *desclcid, elemtype te);
void write_gid_elemtensor(FILE *out, strastre q, long lcid, long sid, const char *desclcid);
void write_gid_elem_type_tensor(FILE *out, strastre q, long lcid, long sid, const char *desclcid, elemtype te);

void write_gid_elem_set(FILE *out, const std::set<elemtype> &reqet, const char *header, long idn1, long ide1);

void export_femcad(FILE *out);
void write_nodes(FILE *out);
void write_elements(FILE *out);
void write_nodes_prep(FILE *out, siftop *st);
void write_elements_prep(FILE *out, siftop *st);
void write_displ(FILE *out, long lcid, const char *deslcid);
void write_nodscalar(FILE *out, strastre scal, long lcid, long dir, const char *desclcid);
void write_nodscalar(FILE *out,double *val, const char *descr, const char *desclcid);
void write_nodscalar(FILE *out, long dir, const char *desclcid);
void write_nforces(FILE *out, long lcid, const char *desclcid, const char *veclabel, double *ifor, bool print_react=true);
void write_elemscalar(FILE *out, strastre scal, long lcid, long dir, const char *desclcid);
void write_elemscalar(FILE *out,double *val, const char *descr, const char *desclcid);
void write_deflection (FILE *out,long lcid,long dir, const char *desclcid);

void print_displacements (FILE *out,long lcid);
void print_multipliers (FILE *out);

void print_strains_nodes (FILE *out,long lcid);
void print_strains_udp (FILE *out,long eid);
//void print_strains_intp (FILE *out,long eid);
//void print_strains (FILE *out,long lcid);

void print_stresses_nodes (FILE *out,long lcid);
void print_stresses_udp (FILE *out,long eid);
//void print_stresses_intp (FILE *out,long lcid);
//void print_stresses (FILE *out,long lcid);

void print_strains_old (FILE *out,long lcid);

void print_stresses_old (FILE *out,long lcid);

void print_other (FILE *out,long lcid);

void print_intforces (FILE *out, double *fi);

void print_forces (FILE *out,double *fi);

void print_reactions (FILE *out,long lcid);

void print_eigenvalues (double *w);

void print_eigenvectors ();

void print_eigenvect_martin (FILE *out);

//void print_edges(FILE *out);

/* termitovo */
void print_valel (gtopology *gt, probdesc *mp, const char *file, char *caption, double *valel,char flag);
//void print_valnod (gtopology *gt, probdesc *mp, mechtop *mt, const char *file, char *caption, double *valnod,char flag);
void print_default_dx (gtopology *gt,probdesc *mp,mechtop *mt,mechmat *mm,long lcid,const char *file);
//void print_default_2_dx (gtopology *gt,probdesc *mp,mechtop *mt,long lcid,const char *file);
/* termitovo */

/// prints header, points and used elements to the VTK file
void print_default_vtk(FILE *out, long lcid, long istep, double lambda);

///prints data to the VTK file
void write_vtk_unkn(FILE *out, long lcid);

void print_vtk_header(FILE *out, long istep, double time);

void aux_mech_nonlin_print (FILE *aux,double *r,double l);
void aux_mech_time_print (FILE *aux,double *r,double l);

#endif
