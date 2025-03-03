#ifndef PGALIAS_H
#define PGALIAS_H
#include "../GEFEL/kwdset.h"

//  aliases for type of domain decomposition method
enum domdectype {no_par_solver=0, schurcompldd=1,fetidd=5,dpfetidd=10,parconjuggrad=40,boundparconjuggrad=45,layered_plate=50};
const enumstr domdectypestr[7] = {{"no_par_solver",0},{"schurcompldd",1}, {"fetidd",5}, {"dpfetidd",10}, {"parconjuggrad",40},{"boundparconjuggrad",45},  {"layered_plate",50}};
const kwdset domdectype_kwdset(sizeof(domdectypestr)/sizeof(*domdectypestr), domdectypestr);

//  aliases for type of solver of coarse problem
/* enum redsystsolver {master_sol=1,paral_sol=2};  */
/* const enumstr redsystsolverstr[2] = {{"master_sol",1}, {"paral_sol",2}};  */
/* const kwdset redsystsolver_kwdset(sizeof(redsystsolverstr)/sizeof(*redsystsolverstr), redsystsolverstr); */

//  aliases for type of preconditioner of FETI method
enum fetiprecond {nofetiprecond=0,lumped=1,dirichlet=2};
const enumstr fetiprecondstr[3] = {{"nofetiprecond",0},{"lumped",1}, {"dirichlet",2}};
const kwdset fetiprecond_kwdset(sizeof(fetiprecondstr)/sizeof(*fetiprecondstr),fetiprecondstr);

// aliases for condensation of fixing node in FETI-DP method
enum condfixing {nocondconer=0,automatic=1,userdef=2};
const enumstr condfixingstr[3] = {{"nocondfixing",0},{"automatic",1}, {"userdef",2}};
const kwdset condfixing_kwdset(sizeof(condfixingstr)/sizeof(*condfixingstr),condfixingstr);

// aliases for method of automatic and user defined condensation of fixing node in FETI-DP method
enum methodcond {nomethod=0,curvecond=1,surfacecond=2,cursurf=3};
const enumstr methodcondstr[4] = {{"nomethod",0},{"curvecond",1},{"surfacecond",2},{"cursurf",3}};
const kwdset methodcond_kwdset(sizeof(methodcondstr)/sizeof(*methodcondstr),methodcondstr);

// aliases for type automatic and user defined condensation of fixing node in FETI-DP method
enum typecondfixing {notype=-1,nocondfixing=0,centroid_fix=1,nth_memb=2,rand_memb=3,n_part_curve=4,all_memb=5,n_mark=6,chose_ring=7,chose_max_ring=8,max_tria=9,userposdef=10};
const enumstr typecondfixingstr[] = {{"notype",-1},{"nocondfixing",0},{"centroid_fix",1},{"nth_memb",2},{"rand_memb",3},{"n_part_curve",4},{"all_memb",5},{"n_mark",6},{"chose_ring",7},{"chose_max_ring",8},{"max_tria",9},{"userposdef",10}}; 
const kwdset typecondfixing_kwdset(sizeof(typecondfixingstr)/sizeof(*typecondfixingstr),typecondfixingstr);


//  aliases for type of preconditioner of parallel conjugate gradient method
enum parcgprec {noprec=0,petscilu=1,pardiagprec=2};
const enumstr parcgprecstr[3] = {{"noprec",0},{"petscilu",1},{"pardiagprec",2}};
const kwdset parcgprec_kwdset(sizeof(parcgprecstr)/sizeof(*parcgprecstr),parcgprecstr);

#endif

