#ifndef PARTALIAS_H
#define PARTALIAS_H
#include "kwdset.h"

enum partType{METIS = 1,JOSTLE = 2,CHACO = 3,PARTY = 4,SCOTCH=5,RALPAR=6};
const enumstr partTypestr[]={{"METIS",1},{"JOSTLE",2},{"CHACO",3},{"PARTY",4},{"SCOTCH",5},{"RALPAR",6}};
const kwdset partType_kwdset(sizeof(partTypestr)/sizeof(*partTypestr), partTypestr);


enum processing{all = 3,preprocessing = 1,postprocessing = 2,aggregates=4};
const enumstr processingstr[]={{"preprocessing",1},{"postprocessing",2},{"all",3},{"aggregates",4}};
const kwdset processing_kwdset(sizeof(processingstr)/sizeof(*processingstr), processingstr);

enum partTech{recursive = 1,kway = 2,vkway = 3,mcrekursive=4,mckway=5,wpartrekursive=6,wpartkway=7,wpartvkway=8};
const enumstr partTechstr[]={{"recursive", 1},{"kway",2},{"vkway",3},{"mcrekursive",4},{"mckway",5},{"wpartrekursive",6},{"wpartkway",7},{"wpartvkway",8}};
const kwdset partTech_kwdset(sizeof(partTechstr)/sizeof(*partTechstr), partTechstr);

enum partOption{defaul = 0,user = 1};
const enumstr partOptionstr[]={{"default", 0},{"user",1}};
const kwdset partOption_kwdset(sizeof(partOptionstr)/sizeof(*partOptionstr), partOptionstr);

enum graphWeight{noweight = 0,nodalweight = 1,edgeweight = 2,nodealedgeweight = 3,multiconstrain=4};
const enumstr graphWeightstr[]={{"noweight",0},{"nodalweight",1},{"edgeweight",2},{"nodealedgeweight",3},{"multiconstrain",4}};
const kwdset graphWeight_kwdset(sizeof(graphWeightstr)/sizeof(*graphWeightstr), graphWeightstr);



#endif
