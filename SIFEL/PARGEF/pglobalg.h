#ifndef PGLOBALG_H
#define PGLOBALG_H

#ifndef EXTERN
#define EXTERN extern
#endif

 //  number of processors
 EXTERN int Nproc;

 //  rank of process
 EXTERN int Myrank;

 //  number of subdomain
 EXTERN int Ndom;

 // name of processors
 EXTERN char proc_name[10000];

 //EXTERN char proc_name[MPI_MAX_PROCESSOR_NAME];
 EXTERN int nameLength;
#endif
