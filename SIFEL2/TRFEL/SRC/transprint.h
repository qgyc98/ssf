#ifndef TRANSPRINT_H
#define TRANSPRINT_H
#include <stdio.h>
#include "gtopology.h"
#include "aliast.h"

void print_initt(long istep, const char *mode, long idn1=1, long ide1=1);
void print_stept(long lcid, long istep, double lambda, double *fi);
void print_stept_forced(long lcid, long istep, double lambda, double *fi);
void print_flusht();
void print_closet();

void export_gid_mesht(FILE *out, long idn1, long ide1);
void export_gid_2dmesht(FILE *out, long icut, long idn1, long ide1);
void export_gid_gauss_ptt(FILE *out);
void write_gid_nodest(FILE *out, long idn1);
void write_gid_elementt(FILE *out, long i, long idn1, long ide1);
void write_gid_2delementt(FILE *out, long i, long id1, long nne, long icut, long di, long idn1, long ide1);
void write_gid_unkn(FILE *out, long lcid, const char *desclcid);
void write_gid_nodvectort(FILE *out, strastret scal, long lcid, long unkn, const char *desclcid);
void write_gid_nforcest(FILE *out, long lcid, const char *desclcid, double *ifor);
void write_gid_nodscalart(FILE *out, strastret scal, long lcid, long dir, const char *desclcid);
void write_gid_elemscalart(FILE *out, strastret scal, long lcid, long dir, const char *desclcid);
void write_gid_elem_type_scalart(FILE *out, strastret scal, long lcid, long dir, const char *desclcid, elemtypet te);
void write_gid_elemvectort(FILE *out, strastret q, long lcid, long dir, const char *desclcid);
void write_gid_elem_type_vectort(FILE *out, strastret q, long lcid, long sid, const char *desclcid, elemtypet te);

void write_nforcest(FILE *out, long lcid, const char *desclcid, double *ifor);
void write_gid_nforcest(FILE* out, long lcid, const char* desclcid, double* ifor);

void print_default_dx (gtopology *gt,long lcid,char *file);
void print_default_vtk (FILE *fname, long lcid, long istep, double lambda);
void write_vtk_unkn(FILE *fname,long lcid, int flag_material, double *fi);
#endif
