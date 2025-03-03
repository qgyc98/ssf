#ifndef ALIAST_H
#define ALIAST_H
#include "kwdset.h"

//  aliases for type of problem solved
enum problemtypet {stationary_problem=50, nonlinear_stationary_problem=51, 
		   nonstationary_problem=60, nonlinear_nonstationary_problem=61, 
		   discont_nonstat_problem=62, discont_nonlin_nonstat_problem=63, 
		   growing_np_problem=70, growing_np_problem_nonlin=71,
		   trfel_tiles=101,hermes=201};
const enumstr problemtypetstr[] = {{"stationary_problem",50},{"nonlinear_stationary_problem",51},
				   {"nonstationary_problem",60},{"nonlinear_nonstationary_problem",61},
				   {"discont_nonstat_problem",62},{"discont_nonlin_nonstat_problem",63},
				   {"growing_np_problem",70},{"growing_np_problem_nonlin",71},
				   {"trfel_tiles",101}, {"hermes",201}};
const kwdset problemtypet_kwdset(sizeof(problemtypetstr)/sizeof(*problemtypetstr),problemtypetstr);


//  aliases for number of transported media/variables
enum transmattert {nomedium=0,onemedium=1,twomediacoup=10,threemediacoup=30,fourmediacoup=40};
const enumstr transmattertstr[] = {{"nomedium",0},{"onemedium",1},{"twomediacoup",10},
				    {"threemediacoup",30},{"fourmediacoup",40}};
const kwdset transmattert_kwdset(sizeof(transmattertstr)/sizeof(*transmattertstr),transmattertstr);


//  aliases for names of transported media/variables
enum mednamest {heat=1,moisture=2,water=3,heat_moisture=10,moisture_salt=20,
		moisture_salt_crystal=25,moisture_salt_crystal_heat=30};
const enumstr mednameststr[] = {{"heat",1},{"moisture",2},{"water",3},{"heat_moisture",10},{"moisture_salt",20},
				 {"moisture_salt_crystal",25},{"moisture_salt_crystal_heat",30}};
const kwdset mednamest_kwdset(sizeof(mednameststr)/sizeof(*mednameststr),mednameststr);



// aliases for name/type of primary variables/unknowns (MUST be numbered from 1 and one by one) 
enum namevart {trf_temperature=1,trf_rel_humidity=2,trf_water_press=3,trf_gas_press=4, trf_pore_press=5,
               trf_salt_conc=6,trf_hydraulic_head=7,trf_vapor_content=8,trf_moisture=9, trf_press_water_vapor=10,trf_capillary_press=3};
const enumstr namevartstr[] = {{"trf_temperature",1}, {"trf_rel_humidity",2}, {"trf_water_press",3},
                               {"trf_gas_press",4},  {"trf_pore_press",5},   {"trf_salt_conc",6},
                               {"trf_hydraulic_head",7},{"trf_vapor_content",8},{"trf_moisture",9}, {"trf_press_water_vapor",10}, {"trf_press_capillary_vapor",11}};
const kwdset namevart_kwdset(sizeof(namevartstr)/sizeof(*namevartstr),namevartstr);


// aliases for solvers of nonstationary problems
enum nonstatsolver {trapezoidd=1,trapezoidv=2,optim_trapezoidd=5,forwardfindif=11};
const enumstr nonstatsolverstr[] = {{"trapezoidd",1}, {"trapezoidv",2}, {"optim_trapezoidd",5},{"forwardfindif",11}};
const kwdset nonstatsolver_kwdset(sizeof(nonstatsolverstr)/sizeof(*nonstatsolverstr), nonstatsolverstr);


//  aliases for nonlinear solvers
enum nlsolvertype {newtont=1};
const enumstr nlsolvertypestr[] = {{"newtont",1}};
const kwdset nlsolvertype_kwdset(sizeof(nlsolvertypestr)/sizeof(*nlsolvertypestr), nlsolvertypestr);



enum gravityaccelerationt {gr_yes=1,gr_no=0};
const enumstr gravityaccelerationtstr[]={{"gr_yes",1},{"gr_no",0}};
const kwdset gravityaccelerationt_kwdset(sizeof(gravityaccelerationtstr)/sizeof(*gravityaccelerationtstr),gravityaccelerationtstr);


enum transpsolver {fullnewtont=1,modnewtont=2};
const enumstr transpsolverstr[]={{"fullnewtont",1},{"modnewtont",2}};
const kwdset transpsolver_kwdset(sizeof(transpsolverstr)/sizeof(*transpsolverstr),transpsolverstr);


enum transpresiduumtype {nonet=0,fluxest=1,lrhst=2};
const enumstr transpresiduumtypestr[]={{"nonet",0}, {"fluxest",1},{"lrhst",2}};
const kwdset transpresiduumtype_kwdset(sizeof(transpresiduumtypestr)/sizeof(*transpresiduumtypestr),transpresiduumtypestr);


enum proptypet {boundarycondt=2, crosssect=4,
		eltypet=6, matelt=7, initcondt=11, comcodnumt=12,
		loadedget=13, loadsource=15, boconsurf=21};


//  aliases for elements
enum elemtypet {barlint=200,barlintax=201,barlint3d=202,barquadt=205,barquadtax=206,trlint=210,trlaxisym=211,
		quadlint=215,quadquadt=216,quadlaxisym=217,quadquadtax=218,
		lineartett=220,linearhext=225,quadratichext=226,linearwedget=228,
		gen2del=300,ifacequadel=351,trquadt=500,trqaxisym=501};//trquadt and trqaxisym not  used
const enumstr elemtypetstr[] = {{"barlint",200},{"barlintax",201},{"barlint3d",202},{"barquadt",205},{"barquadtax",206},
				{"trlint",210},{"trlaxisym",211},{"trquadt",500},{"trqaxisym",501},
				{"quadlint",215},{"quadquadt",216},{"quadlaxisym",217},{"quadquadtax",218},
				{"lineartett",220},{"linearhext",225},{"quadratichext",226},{"linearwedget",228},{"gen2del",300},
                                {"ifacequadel",351}};
const kwdset elemtypet_kwdset(sizeof(elemtypetstr)/sizeof(*elemtypetstr), elemtypetstr);


//  aliases for material models
enum mattypet {isotransmat=100, cernyconcrete=101, nlisotransmat=102,
	       damisotransmat=103, homomat=110, lincoupledmat=120, interfacem=131,
	       bazantpedersen=150, pedersen=151, milly=154, kunzel=155, grunewald=156, devries=157, kunzel2=158, simplediscmat=159,
	       concreteB=160, baroghelB=161, C60baroghelB=165, C30baroghelB=166,
		   o30bazantB=167, C60bazantB=168, C30bazantB=169, glasgow=170, carb1mat=171, moistheat=180,
	       richardsmat=181,
	       salt1mat=200, salt2mat=201, salt3mat=202, salt4mat=203, radiationmater=251,
	       soilmat1=300, discontisotrmat=401, tdisotransmat=421,
	       cementhydrmat=500, sejtkr=650, consolawf1=600, consolwf1=601, consolwf2=700, consolawf2=770, consolhawf3=800, consolhwf2=900};
const enumstr mattypetstr[] = {{"isotransmat",100},{"cernyconcrete",101},{"nlisotransmat",102},
			       {"damisotransmat",103}, {"homomat",110},{"lincoupledmat",120},{"interfacem",131},
			       {"bazantpedersen",150},{"pedersen",151}, {"milly",154}, {"kunzel",155}, {"grunewald",156},
			       {"devries",157}, {"kunzel2",158}, {"simplediscmat",159},
			       {"concreteB",160},{"baroghelB",161},{"C60baroghelB",165},{"C30baroghelB",166},
			       {"o30bazantB",167},{"C60bazantB",168},{"C30bazantB",169},{"glasgow",170},{"carb1mat",171}, {"moistheat",180},
			       {"richardsmat",181},
			       {"salt1mat",200},{"salt2mat",201},{"salt3mat",202},{"salt4mat",203},
			       {"radiationmater",251},
			       {"soilmat1",300},{"discontisotrmat",401},{"tdisotransmat",421},
			       {"cementhydrmat",500},{"sejtkr",650},{"consolawf1",600},{"consolwf1",601},{"consolwf2",700},{"consolawf2",770},{"consolhawf3",800},{"consolhwf2",900}};
const kwdset mattypet_kwdset(sizeof(mattypetstr)/sizeof(*mattypetstr), mattypetstr);


//  aliases for isotherms
enum isotypet {noisotherm=0,ithdata=1,grunewaldroot=2,hansen=3};
const enumstr isotypetstr[] = {{"noisotherm",0},{"ithdata",1},{"grunewaldroot",2},{"hansen",3}};
const kwdset isotypet_kwdset(sizeof(isotypetstr)/sizeof(*isotypetstr), isotypetstr);


//  aliases for cross sections
enum crsectypet {nocrosssectiont=0,crsec1dt=1,crsec2dt=2,crsec3dt=3};
const enumstr crsectypetstr[] = {{"nocrosssectiont",0},{"crsec1dt",1},{"crsec2dt",2},{"crsec3dt",3}};
const kwdset crsectypet_kwdset(sizeof(crsectypetstr)/sizeof(*crsectypetstr), crsectypetstr);


//  aliases for boundary conditions
enum bocontypet {no_bc=0, presc_comp=1, presc_flux=2, det_climcond=3, gen_climcond=5, presc_trmiss=30, presc_trmiss_spec=31,presc_trmiss_isoheat=32,presc_trmiss_computed=50,presc_trmiss_salt=51,presc_trmiss_radiation=90, presc_trmiss_free_surface=40,presc_combined=100};
const enumstr bocontypetstr[] = {{"no_bc",0}, {"presc_comp",1}, {"presc_flux",2}, {"det_climcond",3}, {"gen_climcond",5}, {"presc_trmiss",30}, {"presc_trmiss_spec",31},{"presc_trmiss_isoheat",32},{"presc_trmiss_computed",50},{"presc_trmiss_salt",51},{"presc_trmiss_radiation",90},{"presc_trmiss_free_surface",40},{"presc_combined",100}};
const kwdset bocontype_kwdset(sizeof(bocontypetstr)/sizeof(*bocontypetstr), bocontypetstr);


//  aliases for position on elements, where values are computed
enum elempositiont {nowheret=0,intptst=1,enodest=2,userdefinedt=3};
const enumstr elempositiontstr[] = {{"nowheret",0},{"intptst",1},{"enodest",2},{"userdefinedt",3}};
const kwdset elempositiont_kwdset(sizeof(elempositiontstr)/sizeof(*elempositiontstr), elempositiontstr);


enum strastret {grad=0,flux=1,othert=2,eqothert=3,pgrad=4,pflux=5};
const enumstr strastretstr[] = {{"grad",0}, {"flux",1}, {"othert",2}, {"eqothert",3}, {"pgrad",4}, {"pflux",5}};
const kwdset strastret_kwdset(sizeof(strastretstr)/sizeof(*strastretstr),strastretstr);

//  aliases for format of graphics output
enum graphftt {grftt_no=0, grftt_open_dx=1, grftt_femcad=2, grftt_gid=3, grftt_gid_sep=4, grftt_vtk=5,
               grftt_gid_vtk=6, grftt_gidsep_vtk=7};
const enumstr graphfttstr[] = {{"grftt_no",0}, {"grftt_open_dx",1}, {"grftt_femcad",2}, 
			       {"grftt_gid",3}, {"grftt_gid_sep",4}, {"grftt_vtk",5}, {"grftt_gid_vtk",6}, 
                               {"grftt_gidsep_vtk",7}};
const kwdset graphftt_kwdset(sizeof(graphfttstr)/sizeof(*graphfttstr), graphfttstr);


//  aliases for type of printed unknown
enum prunkt {pr_unknowns=1,pr_gradients=2,pr_fluxes=3,pr_timet=5, pr_stepidt=6,pr_apploadt=7,pr_othert=8,pr_eqothert=9,pr_surffluxes=10};
const enumstr prunktstr[] = {{"pr_unknowns",1}, {"pr_gradients",2}, {"pr_fluxes",3}, {"pr_timet",5},
			      {"pr_stepidt",6}, {"pr_apploadt",7}, {"pr_othert",8}, {"pr_eqothert",9}, {"pr_surffluxes",10}};
const kwdset prunkt_kwdset(sizeof(prunktstr)/sizeof(*prunktstr), prunktstr);


//  aliases for time printing
enum timetypeprint {secondst=1,minutest=2,hourst=3,dayst=4};
const enumstr timetypeprintstr[] = {{"secondst",1}, {"minutest",2}, {"hourst",3}, {"dayst",4}};
const kwdset timetypeprint_kwdset(sizeof(timetypeprintstr)/sizeof(*timetypeprintstr), timetypeprintstr);

//  aliases for type of heat source
enum sourcetype {matfunction=1,concrete_heat=2,cement_hydration=5,seebeck=6};
const enumstr sourcetypestr[] = {{"matfunction",1}, {"concrete_heat",2}, {"cement_hydration",5},{"seebeck",6}};
const kwdset sourcetype_kwdset(sizeof(sourcetypestr)/sizeof(*sourcetypestr),sourcetypestr);

//  aliases for models of water (one phase) flow in deforming porous medium (soils) consol_awf1 = consolidation
enum waterflowtype {lewis_and_schrefler=1,gardner_exponential=2,potts_log_linear=3,van_genuchten=4,consol_camclay=5,consol_camclay_mech=6,lewis_and_schrefler_mefel=7,lewis_and_schrefler_mefel_table=8};
const enumstr waterflowtypestr[] = {{"lewis_and_schrefler",1}, {"gardner_exponential",2}, {"potts_log_linear",3}, {"van_genuchten",4}, {"kuklik_camclay",5}, {"kuklik_camclay_mech",6},{"lewis_and_schrefler_mefel",7},{"lewis_and_schrefler_mefel",8}};
const kwdset waterflowtype_kwdset(sizeof(waterflowtypestr)/sizeof(*waterflowtypestr),waterflowtypestr);

//  aliases for models of air and water (two phase) flow in deforming porous medium (soils) consol_awf2 = consolidation
enum airwaterflowtype {lewis_and_schrefler2=1,van_genuchten2=4,lewis_and_schrefler2_mefel=7};
const enumstr airwaterflowtypestr[] = {{"lewis_and_schrefler2",1}, {"van_genuchten2",4},{"lewis_and_schrefler2_mefel",7}};
const kwdset airwaterflowtype_kwdset(sizeof(airwaterflowtypestr)/sizeof(*airwaterflowtypestr),airwaterflowtypestr);

//  aliases for models of heat, air and water (three phase) flow in deforming porous medium (soils) consol_hawf3 = consolidation
enum heatairwaterflowtype {lewis_and_schrefler3=1,lewis_and_schrefler3_2=2,artificial3=3,van_genuchten3=4,lewis_and_schrefler3_mefel=7};
const enumstr heatairwaterflowtypestr[] = {{"lewis_and_schrefler3",1},{"lewis_and_schrefler3_2",2},{"artificial3",3}, {"van_genuchten3",4},{"lewis_and_schrefler3_mefel",7}};
const kwdset heatairwaterflowtype_kwdset(sizeof(heatairwaterflowtypestr)/sizeof(*heatairwaterflowtypestr),heatairwaterflowtypestr);


//  aliases for models of heat and water (two phase) flow in deforming porous medium (soils) consol_hwf2 = consolidation
enum heatwaterflowtype {lewis_and_schrefler2hw=1,lewis_and_schrefler2hw_2=2,artificial2=3,van_genuchten2hw=4,lewis_and_schrefler2hw_mefel=7};
const enumstr heatwaterflowtypestr[] = {{"lewis_and_schrefler2hw",1},{"lewis_and_schrefler2hw_2",2},{"artificial2",3},{"van_genuchten2hw",4},{"lewis_and_schrefler2hw_mefel",7}};
const kwdset heatwaterflowtype_kwdset(sizeof(heatwaterflowtypestr)/sizeof(*heatwaterflowtypestr),heatwaterflowtypestr);

//  aliases for capacity/reaction computation
enum crcomputation {capacity=1,reaction=2};
const enumstr crcomputationstr[] = {{"capacity",1},{"reaction",2}};
const kwdset crcomputation_kwdset(sizeof(crcomputationstr)/sizeof(*crcomputationstr), crcomputationstr);

enum sourceloc {nod=1,elem=2};

//  aliases for models of retention curve type
enum sr_type {lewis_and_schrefler_sr=1,gardner_exponential_sr=2,potts_log_linear_sr=3,van_genuchten_sr=4,van_genuchten2_sr=5,mefel_sr=7,table_sr=8,bazant_sr=9,baroghel_sr=10,masin_sr=11,febex_granit_sr=12};
const enumstr sr_typestr[] = {{"lewis_and_schrefler_sr",1}, {"gardner_exponential_sr",2}, {"potts_log_linear_sr",3}, {"van_genuchten_sr",4}, 
{"mefel_sr",7},{"table_sr",8},{"bazant_sr",9},{"baroghel_sr",10},{"masin_sr",11},{"febex_granit_sr",12}};
const kwdset sr_type_kwdset(sizeof(sr_typestr)/sizeof(*sr_typestr),sr_typestr);

//  aliases for models of effective stress factor (parameter) type
enum xi_type {biot_xi=1,masin_xi=2,biot_reduced_xi=3,biot_masin_xi=4,masin_hypopl_xi=5};
const enumstr xi_typestr[] = {{"biot_xi",1}, {"masin_xi",2},{"biot_reduced_xi",3},{"biot_masin_xi",4},{"masin_hypopl_xi",5}};
const kwdset xi_type_kwdset(sizeof(xi_typestr)/sizeof(*xi_typestr),xi_typestr);


#endif

