#ifndef __IO_H 
#define __IO_H

class IO {

public: 
	long nn;
	long ne;
	long *nadjnodnod;		// pocet sousedicich bodu s bodem
	long **adjnodnod;		// cisla sousedicich bodu s bodem 
	long *nadjelnod;
	long **adjelnod;
	long *nadjelel;
	long **adjelel;
	long * nnod;
	long ** nodes;		// element = seznam cisel uzlu
	double **xyz;			// bod = souradnice x,y,z
	IO (const char* fileName);
	~IO();
};

#endif

