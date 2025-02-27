#ifndef CONVERTORALIAS_H
#define CONVERTORALIAS_H
#include "kwdset.h"

enum convType{T3D2T3D = 1,T3D2SIFEL = 2,T3D2GID = 3,SIFEL2GID};
const enumstr convTypestr[]={{"T3D2T3D",1},{"T3D2SIFEL",2},{"T3D2GID",3},{"SIFEL2GID",4}};
const kwdset convType_kwdset(sizeof(convTypestr)/sizeof(*convTypestr), convTypestr);


enum convprocessing{seq = 1,paral = 2};
const enumstr convprocessingstr[]={{"seq",1},{"paral",2}};
const kwdset convprocessing_kwdset(sizeof(convprocessingstr)/sizeof(*convprocessingstr), convprocessingstr);




#endif
