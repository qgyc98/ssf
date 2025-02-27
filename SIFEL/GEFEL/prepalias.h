#ifndef PREPALIAS_H
#define PREPALIAS_H

#include "kwdset.h"

// keywords for beginnings of preprocessor sections
enum  bsec {begsec_files=0,    begsec_probdesc=1, begsec_loadcase=2, 
            begsec_nodvertpr=3,begsec_nodedgpr=4, begsec_nodsurfpr=5, begsec_nodvolpr=6,
            begsec_eledgpr=7,  begsec_elsurfpr=8, begsec_elvolpr=9,
            begsec_outdrv=10,  begsec_part=12,    begsec_convert=13,  begsec_mater=14,
            begsec_crsec=15, begsec_elvertpr=16};
const enumstr bsec_str[] = {{"begsec_files",0},     {"begsec_probdesc",1}, {"begsec_loadcase",2},
                            {"begsec_nodvertpr",3}, {"begsec_nodedgpr",4}, {"begsec_nodsurfpr", 5}, 
                            {"begsec_nodvolpr",6},  {"begsec_eledgpr",7},  {"begsec_elsurfpr",8}, 
                            {"begsec_elvolpr",9},   {"begsec_outdrv",10},  {"begsec_gfunct", 11}, 
                            {"begsec_part",12},     {"begsec_convert",13}, {"begsec_mater",14},    
                            {"begsec_crsec",15},    {"begsec_elvertpr",16}};
const kwdset bsec_kwdset(sizeof(bsec_str)/sizeof(*bsec_str), bsec_str);

// keywords for ends of preprocessor sections
enum  esec {endsec_files=0,    endsec_probdesc=1, endsec_loadcase=2, 
            endsec_nodvertpr=3,endsec_nodedgpr=4, endsec_nodsurfpr=5, endsec_nodvolpr=6,
            endsec_eledgpr=7,  endsec_elsurfpr=8, endsec_elvolpr=9,   endsec_outdrv=10,  
            endsec_part=12,    endsec_convert=13, endsec_mater=14,    endsec_crsec=15,
            endsec_elvertpr=16};
const enumstr esec_str[] = {{"endsec_files",0},     {"endsec_probdesc",1}, {"endsec_loadcase",2},
                            {"endsec_nodvertpr",3}, {"endsec_nodedgpr",4}, {"endsec_nodsurfpr", 5}, 
                            {"endsec_nodvolpr",6},  {"endsec_eledgpr",7},  {"endsec_elsurfpr",8}, 
                            {"endsec_elvolpr",9},   {"endsec_outdrv",10},  {"endsec_gfunct",11},  
                            {"endsec_part",12},     {"endsec_convert",13}, {"endsec_mater",14},   
                            {"endsec_crsec",15},    {"endsec_elvertpr",16}};
const kwdset esec_kwdset(sizeof(esec_str)/sizeof(*esec_str), esec_str);

enum  prepcondtype {cond=1, file=2};
const enumstr prepcondtype_str[] = {{"cond",1}, {"file",2}};
const kwdset prepcondtype_kwdset(sizeof(prepcondtype_str)/sizeof(*prepcondtype_str), prepcondtype_str);

#endif
