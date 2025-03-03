#ifndef ALIASC_H
#define ALIASC_H
#include "kwdset.h"

enum problemtypec {fully_coupled_mech_trans=1,fully_coupled_material=2,
		   par_coupl_mech_trans=10,growing_par_coupl_mech_trans=70};
const enumstr problemtypecstr[4] = {{"fully_coupled_mech_trans",1},{"fully_coupled_material",2},
				    {"par_coupl_mech_trans",10},{"growing_par_coupl_mech_trans",70}};
const kwdset problemtypec_kwdset(sizeof(problemtypecstr)/sizeof(*problemtypecstr),problemtypecstr);


enum babuskabrezzi {notdefbb=0,lin_lin=1,quad_lin=2,quad_quad=3};
const enumstr babuskabrezzistr[4] = {{"notdefbb",0},{"lin_lin",1},{"quad_lin",2},{"quad_quad",3}};
const kwdset babuskabrezzi_kwdset(sizeof(babuskabrezzistr)/sizeof(*babuskabrezzistr),babuskabrezzistr);

enum coupcleanmatrices {clean_no=0, clean_yes=1};
const enumstr coupcleanmatricesstr[2] = {{"clean_no",0},{"clean_yes",1}};
const kwdset coupcleanmatrices_kwdset(sizeof(coupcleanmatricesstr)/sizeof(*coupcleanmatricesstr),coupcleanmatricesstr);

enum coupsolver {linearc=0, fullnewtonc=1, modnewtonc=2};
const enumstr coupsolverstr[3] = {{"linearc",0},{"fullnewtonc",1},{"modnewtonc",2}};
const kwdset coupsolver_kwdset(sizeof(coupsolverstr)/sizeof(*coupsolverstr),coupsolverstr);

enum residuumtype {fluxesc=1,lrhsc=2};
const enumstr residuumtypestr[2]={{"fluxesc",1},{"lhrsc",2}};
const kwdset residuumtype_kwdset(sizeof(residuumtypestr)/sizeof(*residuumtypestr),residuumtypestr);

enum transmatterc {mech_onemedium=1,mech_twomedia=2,mech_threemedia=3,mech_fourmedia=4};
const enumstr transmattercstr[4] = {{"mech_onemedium",1},{"mech_twomedia",2},{"mech_threemedia",3},{"mech_fourmedia",4}};
const kwdset transmatterc_kwdset(sizeof(transmattercstr)/sizeof(*transmattercstr),transmattercstr);

enum mednamesc {mech_heat=1,mech_moisture=2,mech_water=3,mech_air_water=4,mech_heat_water=5,mech_heat_moisture=10,mech_moisture_salt=15,mech_heat_moisture_salt=20};
const enumstr mednamescstr[8] = {{"mech_heat",1},{"mech_moisture",2},{"mech_water",3},{"mech_air_water",4},{"mech_heat_water",5},{"mech_heat_moisture",10},{"mech_moisture_salt",15},{"mech_heat_moisture_salt",20}};
const kwdset mednamesc_kwdset(sizeof(mednamescstr)/sizeof(*mednamescstr),mednamescstr);


enum gravityaccelerationc {grc_yes=1,grc_no=0};
const enumstr gravityaccelerationcstr[]={{"grc_yes",1},{"grc_no",0}};
const kwdset gravityaccelerationc_kwdset(sizeof(gravityaccelerationcstr)/sizeof(*gravityaccelerationcstr),gravityaccelerationcstr);

enum nonlinsolvertypec {arclc=1,newtonc=2,newtondform=3};
const enumstr nonlinsolvertypecstr[3] = {{"arclc",1},{"newtonc",2},{"newtondform",3}};
const kwdset nonlinsolvertypec_kwdset(sizeof(nonlinsolvertypecstr)/sizeof(*nonlinsolvertypecstr),nonlinsolvertypecstr);


enum displacementnormc {alldofsc=1,seldofsc=2,selnodesc=3,nodesdistincrc=4};


enum elemtypec {coupbar=10,coupquad=20,couptria=30,coupaxiquad=40,axisymfc=45,coupaxitria=50,couphex=60,couptet=70};
const enumstr elemtypecstr[8] = {{"coupbar",10},{"coupquad",20},{"couptria",30},
				 {"coupaxiquad",40},{"axisymfc",45},{"coupaxitria",50},{"couphex",60},{"couptet",70}};
const kwdset elemtypec_kwdset(sizeof(elemtypecstr)/sizeof(*elemtypecstr),elemtypecstr);


enum mattypec {isotransmatc=100,sejtkrc=650,concreteBc=160,baroghelBc=161,C60baroghelBc=165,
               C30baroghelBc=166,o30bazantBc=167,C60bazantBc=168,glasgowc=170,glascoup=171,
	       bentonite=400,
	       consolawf1c=600,consolwf1c=601,consolwf2c=700,consolawf2c=770,consolhawf3c=800,consolhwf2c=900};
const enumstr mattypecstr[17] = {{"isotransmatc",100},{"sejtkrc",650},{"concreteBc",160},{"baroghelBc",161},{"C60baroghelBc",165},
				 {"C30baroghelBc",166},{"o30bazantBc",167},{"C60bazantBc",168},{"glasgowc",170},{"glascoup",171},
				 {"bentonite",400},
				 {"consolawf1c",600}, {"consolwf1c",601}, {"consolwf2c",700}, {"consolawf2c",770},{"consolhawf3c",800},{"consolhwf2c",900}};
const kwdset mattypec_kwdset(sizeof(mattypecstr)/sizeof(*mattypecstr),mattypecstr);


// aliases for keyword processing switch
enum pkwd_swc {nokwdc=0, kwd_pdc=1, kwd_outdc=2, kwd_pd_outdc=3, kwd_fullc=4};

//  aliases for models of water (one phase) flow in deforming porous medium (soils) consol_awf1 = consolidation
enum waterflowmechtype {lewis_and_schrefler_coup=1,van_genuchten_coup=4,consol_camclay_coup=5,consol_camclay_mech_coup=6,lewis_and_schrefler_mefel_coup=7};
const enumstr waterflowmechtypestr[] = {{"lewis_and_schrefler_coup",1},{"van_genuchten_coup",4},{"consol_camclay_coup",5},{"consol_camclay_mech_coup",6},{"lewis_and_schrefler_mefel_coup",7}};
const kwdset waterflowmechtype_kwdset(sizeof(waterflowmechtypestr)/sizeof(*waterflowmechtypestr),waterflowmechtypestr);

//  aliases for models of air and water (two phase) flow in deforming porous medium (soils) consol_awf2 = consolidation
enum airwaterflowmechtype {lewis_and_schrefler2_coup=1,van_genuchten2_coup=4,lewis_and_schrefler2_mefel_coup=7};
const enumstr airwaterflowmechtypestr[] = {{"lewis_and_schrefler2",1},{"van_genuchten2_coup",4},{"lewis_and_schrefler2_mefel",7}};
const kwdset airwaterflowmechtype_kwdset(sizeof(airwaterflowmechtypestr)/sizeof(*airwaterflowmechtypestr),airwaterflowmechtypestr);

//  aliases for models of heat and water (two phase) flow in deforming porous medium (soils) consol_hwf2 = consolidation
enum heatwaterflowmechtype {lewis_and_schrefler2hw_coup=1,lewis_and_schrefler2hw_mefel_coup=7};
const enumstr heatwaterflowmechtypestr[] = {{"lewis_and_schrefler2hw",1},{"lewis_and_schrefler2hw_mefel",7}};
const kwdset heatwaterflowmechtype_kwdset(sizeof(heatwaterflowmechtypestr)/sizeof(*heatwaterflowmechtypestr),heatwaterflowmechtypestr);

//  aliases for models of heat, air and water (three phase) flow in deforming porous medium (soils) consol_hawf3 = consolidation
enum heatairwaterflowmechtype {lewis_and_schrefler3_coup=1,lewis_and_schrefler3_2_coup=2,lewis_and_schrefler3_mefel_coup=7};
const enumstr heatairwaterflowmechtypestr[] = {{"lewis_and_schrefler3",1},{"lewis_and_schrefler3_2",2},{"lewis_and_schrefler3_mefel",7}};
const kwdset heatairwaterflowmechtype_kwdset(sizeof(heatairwaterflowmechtypestr)/sizeof(*heatairwaterflowmechtypestr),heatairwaterflowmechtypestr);

// aliases for the various data passing approaches between modules
// pass_by_closest_ip=1 -> data are copied to nodes from the closest int. point, requires the same number of elements in both modules - minimum additional memory requirements, fastest way, may be inaccurate for certain combination of elements
// pass_by_nodes_comp=2 -> data are calculated at nodes with the help of state varibales from the closest int. point, requires the first n elements in both modules to be the same in shape and same ordering of the first m nodes,
//                         where n is the minimum number of elements in modules, m is the minimum number of nodes on elements - minimum additional memory requirements, slower than type 1, better accuracy than type 1
// pass_by_aux_ip = 3 -> data are calculated in auxiliary int. points, the number of nodes and elements in meshes of particular modules can be independent - may leed to large memory requirements, slower than type 1, exact.
// pass_by_copy_ip = 4 -> data are copied directly from integration points (the SAME int. points must be used for the first Tt->ne elements on both meshes) - exact, fast, minimum memory requirements.
enum datapasstype {pass_by_closest_ip=1, pass_by_nodes_comp=2, pass_by_aux_ip=3, pass_by_copy_ip=4};
const enumstr datapasstypestr[] = {{"pass_by_closest_ip", 1}, {"pass_by_nodes_comp", 2}, {"pass_by_aux_ip", 3}, {"pass_by_copy_ip", 4}};
const kwdset datapasstype_kwdset(sizeof(datapasstypestr)/sizeof(*datapasstypestr),datapasstypestr);
#endif


