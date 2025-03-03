#ifndef ALIAS_H
#define ALIAS_H
#include "kwdset.h"

// aliases for type of solved problem
enum problemtype {linear_statics=1, eigen_dynamics=2, forced_dynamics=3,
                  linear_stability=5,mat_nonlinear_statics=10,geom_nonlinear_statics=11,
                  mech_timedependent_prob=15,growing_mech_structure=17,gmech_timedependent_prob=18,
                  earth_pressure=20,layered_linear_statics=30,
                  lin_floating_subdomain=35,nonlin_floating_subdomain=36,
                  hemivar_inequal=37,
                  nonlinear_dynamics=40,load_balancing=51};
const enumstr problemtypestr[] = {{"linear_statics",1}, {"eigen_dynamics",2}, {"forced_dynamics",3},
                                  {"linear_stability",5}, {"mat_nonlinear_statics",10}, {"geom_nonlinear_statics",11},
                                  {"mech_timedependent_prob",15}, {"growing_mech_structure",17}, {"gmech_timedependent_prob",18},
                                  {"earth_pressure",20}, {"var_stiff_method",25}, 
                                  {"layered_linear_statics",30}, 
                                  {"lin_floating_subdomain",35}, {"nonlin_floating_subdomain",36}, 
                                  {"hemivar_inequal",37}, {"nonlinear_dynamics",40}, {"load_balancing",51}};
const kwdset problemtype_kwdset(sizeof(problemtypestr)/sizeof(*problemtypestr), problemtypestr);



// aliases for forced dynamics solvers
enum forcedsolver {newmark=1,findiff=2,explfindiff=3,modal_analysis=51};
const enumstr forcedsolverstr[] = {{"newmark",1}, {"findiff",2}, {"explfindiff",3}, {"modal_analysis",51}};
const kwdset forcedsolver_kwdset(sizeof(forcedsolverstr)/sizeof(*forcedsolverstr), forcedsolverstr);



//aliases for nonlinear statics solvers
enum nonlinsolvertype {arcl=1,newton=2,arclrv1=3,newtonrv1=4,newtonrestart=5,displctrl=6,displctrlrv=7,arclrv=8,adaptram=10,arclg=20,newtong=30,dissip_incr=40};
const enumstr nonlinsolvertypestr[] = {{"arcl",1}, {"newton",2}, {"arclrv1",3}, {"newtonrv1",4}, {"newtonrestart",5}, 
                                       {"displctrl",6}, {"displctrlrv",7}, {"arclrv",8}, {"adaptram",10}, {"arclg",20}, {"newtong", 30}, {"dissip_incr", 40}};
const kwdset nonlinsolvertype_kwdset(sizeof(nonlinsolvertypestr)/sizeof(*nonlinsolvertypestr), nonlinsolvertypestr);



//aliases for nonlinear time-dependent solvers
enum nonlintimesolvertype {fullnewton=1,modnewton=2};
const enumstr nonlintimesolvertypestr[] = {{"fullnewton",1},{"modnewton",2}};
const kwdset  nonlintimesolvertype_kwdset(sizeof(nonlintimesolvertypestr)/sizeof(*nonlintimesolvertypestr), nonlintimesolvertypestr);



// aliases for earth pressure solvers
enum epsolvertype {gep_sol=1, gepvarsup_sol=2, epplast_sol=3};
const enumstr epsolvertypestr[] = {{"gep_sol",1}, {"gepvarsup_sol",2}, {"epplast_sol",3}};
const kwdset epsolvertypes_kwdset(sizeof(epsolvertypestr)/sizeof(*epsolvertypestr), epsolvertypestr);



// aliases for type of displacement norm computation
enum displacementnorm {alldofs=1,seldofs=2,seldofscoord=3, selmstr=4, selecnodes=6,nodesdistincr=8};
const enumstr displacementnormstr[] = {{"alldofs",1}, {"seldofs",2}, {"seldofscoord",3}, {"selmstr", 4}, {"selecnodes",6}, {"nodesdistincr",8}};
const kwdset displacementnorm_kwdset(sizeof(displacementnormstr)/sizeof(*displacementnormstr), displacementnormstr);



// aliases for types for stress return algorithms
enum stressretalgtype {nostressretalg=0,cp=1,gsra=10,rkfsra=11};
const enumstr stressretalgtypestr[] = {{"nostressretalg",0}, {"cp",1}, {"gsra",10}, {"rkfsra",11}};
const kwdset stressretalgtype_kwdset(sizeof(stressretalgtypestr)/sizeof(*stressretalgtypestr), stressretalgtypestr);


// aliases for types for stress return algorithms of Runge-Kutta type
enum rktype {nork=0, fwdeulert=1, heunt=2, rkf23t=3, rkf23bst=4, rkf34t=5, rkf45t=6, fwdeulert_dirstiff=11, heunt_dirstiff=12, rkf23t_dirstiff=13, rkf23bst_dirstiff=14, rkf34t_dirstiff=15, rkf45t_dirstiff=16};
const enumstr rktypestr[] = {{"nork",0}, {"fwdeulert",1}, {"heunt",2}, {"rkf23t",3}, {"rkf23bst",4}, {"rkf34t",5}, {"rkf45t",6}, {"fwdeulerfsub",7}, {"fwdeulert_dirstiff",11}, {"heunt_dirstiff",12}, {"rkf23t_dirstiff",13}, {"rkf23bst_dirstiff",14}, {"rkf34t_dirstiff",15}, {"rkf45t_dirstiff",16}, {"fwdeulerfsub_dirstiff",17}};
const kwdset rktype_kwdset(sizeof(rktypestr)/sizeof(*rktypestr), rktypestr);


// aliases for types of stiffness matrix computation
enum stiffmatrix {initial_stiff=1,tangent_stiff=2,secant_stiff=3,incr_tangent_stiff=5,ijth_tangent_stiff=6};
const enumstr stiffmatrixstr[] = {{"initial_stiff",1}, {"tangent_stiff",2}, {"secant_stiff",3}, {"incr_tangent_stiff",5}, {"ijth_tangent_stiff",6}};
const kwdset stiffmatrix_kwdset(sizeof(stiffmatrixstr)/sizeof(*stiffmatrixstr), stiffmatrixstr);

// aliases for determination of increment of lambda in arc-length
enum detlambda {nodetermination=0, minvalue=1, maxvalue=2, minangle=3, linearizedmeth=4, fullmethod=5};
const enumstr detlambdastr[] = {{"nodetermination",0}, {"minvalue",1}, {"maxvalue",2}, {"minangle",3}, {"linearizedmeth",4}, {"fullmethod",5}};
const kwdset detlambda_kwdset(sizeof(detlambdastr)/sizeof(*detlambdastr), detlambdastr);


// aliases for damping types
enum damping {nodamp=0,massdamp=1,stiffdamp=2,rayleighdamp=3,elemdamp=11};
const enumstr dampingstr[] = {{"nodamp",0}, {"massdamp",1}, {"stiffdamp",2}, {"rayleighdamp",3}, {"elemdamp",11}};
const kwdset damping_kwdset(sizeof(dampingstr)/sizeof(*dampingstr), dampingstr);



// aliases used for properties in preprocessor
enum proptype {ndofs=1, boundarycond=2, loadnodes=3, crosssec=4, localsystem=5,
               eltype=6, matel=7, loadelems=8, dloadnodes=9, dloadelems=10,
               initcond=11, comcodnum=12, loadedge=13, sscomp=14, nodetemp=15, 
               timefun=16,surfaceload=17, elsurfaceload=18, nodtimefunc=19, 
               eltimefunc=20};



// aliases for element types
enum elemtype {bar2d=1,beam2d=2,bar3d=3,beam3d=4,beamg3d=5,barq2d=6,barq3d=7,subsoilbeam=8,beam2dsp=9,
               spring_1 = 10, spring_2 = 11, spring_3 = 12,
               spring_4 = 13, spring_5 = 14, spring_6 = 15,
               planeelementlt=20,planeelementqt=21,planeelementrotlt=22,
               planeelementlq=23,planeelementqq=24,planeelementrotlq=25,
               planeelementsubqt=30,planequadinterface=35,
               cctel=41,dktel=42,dstel=43,q4plateel=45,argyristr=46,quadkirch=47,dkqel=48,
	       subsoilplatetr=50,subsoilplateq=51,
               axisymmlt=60,axisymmqt=61,axisymmlq=63,axisymmqq=64,axisymmcq=65, axisymmlqintface=66,
               shelltrelem=80,shellqelem=81,shelltrmelem=82,
               lineartet=100,quadrtet=101,linearhex=102,quadrhex=103,lineartetrot=104,linearhexrot=105,
	       linearwed=106,quadrwed=107,hexintface=111,
               particleelem=200,tetralatt=250};
const enumstr elemtypestr[] = {{"bar2d",1}, {"beam2d",2}, {"bar3d",3}, {"beam3d",4}, {"beamg3d",5}, 
                               {"barq2d",6},{"barq3d",7}, {"subsoilbeam",8},
                               {"spring_1",10}, {"spring_2",11}, {"spring_3",12},
                               {"spring_4",13}, {"spring_5",14}, {"spring_6",15},
                               {"planeelementlt",20}, {"planeelementqt",21}, {"planeelementrotlt",22},
                               {"planeelementlq",23}, {"planeelementqq",24}, {"planeelementrotlq",25},
                               {"planeelementsubqt",30}, {"planequadinterface",35},
                               {"cctel",41}, {"dktel",42}, {"dstel",43}, {"q4plateel",45}, {"argyristr",46},
                               {"quadkirch",47}, {"dkqel",48}, {"subsoilplatetr",50}, {"subsoilplateq",51},
                               {"axisymmlt",60}, {"axisymmqt",61}, {"axisymmlq",63}, {"axisymmqq",64}, {"axisymmcq",65},
                               {"axisymmlqintface",66},
                               {"shelltrelem",80}, {"shellqelem",81}, {"shelltrmelem",82},
                               {"lineartet",100},  {"quadrtet",101}, {"linearhex",102}, {"quadrhex",103},
                               {"lineartetrot",104}, {"linearhexrot",105},
                               {"linearwed",106}, {"quadrwed",107}, {"hexintface",111},
                               {"particleelem",200}, {"tetralatt",250}};
const kwdset elemtype_kwdset(sizeof(elemtypestr)/sizeof(*elemtypestr), elemtypestr);



// aliases for strain/stress state
enum strastrestate {bar=1,plbeam=2,spacebeam=5,
                    planestress=10,planestrain=11,planecontact=12,platek=15,plates=16,
                    axisymm=20,shell=25,
                    spacestress=30};
const enumstr strastrestatestr[] = {{"bar",1}, {"plbeam",2}, {"spacebeam",5},
                                    {"planestress",10}, {"planestrain",11}, {"planecontact",12}, {"platek",15}, {"plates",16},
                                    {"axisymm",20}, {"shell",25}, {"spacestress",30}};
const kwdset strastrestate_kwdset(sizeof(strastrestatestr)/sizeof(*strastrestatestr), strastrestatestr);



// aliases for  material model classes
enum materialmodel {local=0,nonlocal=1};
const enumstr materialmodelstr[] = {{"local",0}, {"nonlocal",1}};
const kwdset materialmodel_kwdset(sizeof(materialmodelstr)/sizeof(*materialmodelstr), materialmodelstr);


// aliases for material models
enum mattype {elisomat=1,elgmat3d=2,elortomat=3,elgmat2d=4, homomatm=9,
              simplas1d=10,jflow=11,mohrcoulparab=22,mohrcoul=24,
              boermaterial=25,druckerprager=26,doubledrprager=27,druckerprager2=28,modcamclaymat=30, modcamclaycoupmat=31,bbmcoupmat=33,doublestructuremat=34,
              shefpl=40,chenplast=42, hissplasmat=45,
              microplaneM4=50,microsimp=51,microfibro=52,
              simvisplas=70,simvisc=71,isovisc=72,lemaitr=75,
              layerplate=80,
              scaldamage=100,scaldamagecc=101,glasgowdamage=102,anisodamage=104, anisodamagerot=105, ortodamage=106, ortodamagerot=107, fixortodamage=108, ortodamage2=109,
              contmat=150, cebfipcontmat=160, damplifmat=170, plastifmat=180,
              aci=200,cebfip=202,
              nonlocplastmat=310,nonlocdamgmat=320,nonlocalmicroM4=340,
              graphm=400, geoelast=401, varelisomat=405, elisopdmat=406, hypoplastmat=420,  hypoplastusatthermat=421, 
              creepbaz=500,creepb3=502,creepdpl=503,creeprs=504,creepeffym=505,winklerpasternak=550,consolidation=600,
              glasgowmechmat=700,
              rcmatmodnorm=800,
              therisodilat=900,thervolisodilat=901,relaxationeuro=951,
              damage_plasticity=1000, viscoplasticity=1010, viscoelasticity=1020,
              creep_damage=1030, time_switchmat=1040, effective_stress=1050, shrinkagemat=1060, elasttimemat=1070,
              lenjonespot=2000, cusatismat=2010};
const enumstr mattypestr[] = {{"elisomat",1}, {"elgmat3d",2}, {"elortomat",3}, {"elgmat2d",4},{"homomatm",9},{"simplas1d",10}, 
                              {"jflow",11},  {"mohrcoulparab",22}, 
                              {"mohrcoul",24}, {"boermaterial",25}, {"druckerprager",26},  {"doubledrprager",27}, {"druckerprager2",28}, 
                              {"modcamclaymat",30}, {"modcamclaycoupmat",31}, {"bbmcoupmat",33}, {"doublestructuremat",34}, {"shefpl",40}, {"chenplast",42}, {"hissplasmat",45}, 
                              {"microplaneM4",50}, {"microsimp",51}, {"microfibro",52}, 
                              {"simvisc",71}, {"isovisc",72}, {"lemaitr",75}, 
                              {"layerplate",80},
                              {"scaldamage",100}, {"scaldamagecc",101}, {"glasgowdamage",102}, 
                              {"anisodamage", 104 }, {"anisodamagerot", 105 }, {"ortodamage", 106 }, { "ortodamagerot", 107 }, {"fixortodamage", 108}, {"ortodamage2", 109},
                              {"contmat",150}, {"cebfipcontmat",160}, {"damplifmat",170}, {"plastifmat",180},
                              {"aci",200}, {"cebfip",202}, {"nonlocplastmat",310}, 
                              {"nonlocdamgmat",320}, {"nonlocalmicroM4",340}, {"graphm",400}, 
                              {"geoelast",401}, {"varelisomat",405}, {"elisopdmat",406}, {"hypoplastmat",420}, {"hypoplastusatthermat",421}, {"creepbaz",500}, 
                              {"creepb3",502}, {"creepdpl",503}, {"creeprs",504}, {"creepeffym",505},
                              {"winklerpasternak",550}, {"consolidation",600}, {"glasgowmechmat",700}, 
                              {"therisodilat",900},{"thervolisodilat",901},{"relaxationeuro",951},
                              {"damage_plasticity",1000}, {"viscoplasticity",1010}, {"viscoelasticity",1020},
                              {"creep_damage",1030}, {"time_switchmat",1040}, {"effective_stress",1050}, {"shrinkagemat",1060},
                              {"elasttimemat",1070}, {"lenjonespot",2000}, {"cusatismat",2010}};
const kwdset mattype_kwdset(sizeof(mattypestr)/sizeof(*mattypestr), mattypestr);



// aliases for cross-section 
enum crsectype {nocrosssection=0, csbar2d=1,csbeam2d=2,csbeam3d=4,csplanestr=10,cs3dprob=20,csnodal=40,cslayer=50};
const enumstr crsectypestr[] = {{"nocrosssection",0}, {"csbar2d",1}, {"csbeam2d",2}, {"csbeam3d",4}, 
                                {"csplanestr",10}, {"cs3dprob",20}, {"csnodal",40}, {"cslayer",50}};
const kwdset crsectype_kwdset(sizeof(crsectypestr)/sizeof(*crsectypestr), crsectypestr);



// aliases for load of element
enum elloadtype {noelload=0, volume=1, edge=2,surface=3,edge_surface=4,edge_volume=5,surface_volume=6, edge_surface_volume=7, beam_load=8};
const enumstr elloadtypestr[] = {{"noelload",0}, {"volume",1}, {"edge",2}, {"surface",3}, {"edge_surface",4}, 
                                 {"edge_volume",5}, {"surface_volume",6}, {"edge_surface_volume",7}, {"beam_load",8}};
const kwdset elloadtype_kwdset(sizeof(elloadtypestr)/sizeof(*elloadtypestr), elloadtypestr);


// aliases for element load meaning
enum elloadmeaning {noelloadmean=0, load_contmech_glob=1, load_contmech_loc=2, load_pointmech_glob=3, load_pointmech_loc=4, load_thermgrad_y=5, load_thermgrad_z=6};
const enumstr elloadmeaningstr[] = {{"noelloadmean",0}, {"load_contmech_glob",1}, {"load_contmech_loc",2}, {"load_pointmech_glob",3}, {"load_pointmech_loc",4}, 
                                    {"load_thermgrad_y",5}, {"load_thermgrad_z",6}};
const kwdset elloadmeaning_kwdset(sizeof(elloadmeaningstr)/sizeof(*elloadmeaningstr), elloadmeaningstr);


// aliases for point position on element
enum elemposition {nowhere=0,intpts=1,enodes=2,cenodes=3,userdefined=10};
const enumstr elempositionstr[] = {{"nowhere",0}, {"intpts",1}, {"enodes",2}, {"cenodes",3}, {"userdefined",10}};
const kwdset elemposition_kwdset(sizeof(elempositionstr)/sizeof(*elempositionstr), elempositionstr);



// aliases for quantity types
enum strastre {undef_strastre=-1,strain=0,stress=1,other=2,pstrain=3, pstress=4};
const enumstr strastrestr[] = {{"undef_strastre",-1}, {"strain",0}, {"stress",1}, {"other",2}, {"pstrain",3}, {"pstress",4}};
const kwdset strastre_kwdset(sizeof(strastrestr)/sizeof(*strastrestr), strastrestr);



// aliases for initial conditions and initial values
enum inictype {none = 0, inidisp = 1, inistrain = 2, inistress = 4, iniother = 8, inicond = 16, inidisp_x = 32, 
               inidisp_y = 64, inidisp_z = 128, iniderdisp = 256, inidisp_derdisp = 257};
const enumstr inictypestr[] = {{"none",0}, {"inidisp",1}, {"inistrain",2}, 
                               {"inistress",4}, {"iniother",8}, {"inicond",16}, {"inidisp_x", 32}, 
                               {"inidisp_y",64}, {"inidisp_z",128}, {"iniderdisp",256}, {"inidisp_derdisp", 257}};
const kwdset inictype_kwdset(sizeof(inictypestr)/sizeof(*inictypestr), inictypestr);



// aliases for types of equivalent strain computation
enum paramf_type {norstrain = 1,    norenergy = 2,
                  norposstrain = 3, norposenergy = 4,
                  norrankine = 5,   norrankinesmooth = 6,
                  normazar = 7, vonmises= 8};
const enumstr paramf_typestr[] = {{"norstrain",1}, {"norenergy",2}, {"norposstrain",3}, {"norposenergy",4},
                                  {"norrankine",5}, {"norrankinesmooth",6}, {"normazar",7}, {"vonmises",8}};
const kwdset paramf_type_kwdset(sizeof(paramf_typestr)/sizeof(*paramf_typestr), paramf_typestr);



// aliases for types of damage function
enum damfunc_type {simpexp = 1,    mazarsexp = 2};
const enumstr damfunc_typestr[] = {{"simpexp",1}, {"mazarsexp",2}};
const kwdset damfunc_type_kwdset(sizeof(damfunc_typestr)/sizeof(*damfunc_typestr), damfunc_typestr);



// aliases for switching on/off of correction of dissipated energy
enum corr_disip_en {corr_off = 0, corr_on = 1}; 
const enumstr corr_disip_enstr[] = {{"corr_off",0}, {"corr_on",1}}; 
const kwdset corr_disip_en_kwdset(sizeof(corr_disip_enstr)/sizeof(*corr_disip_enstr), corr_disip_enstr);



// aliases for evolution function of damage parmetrs in anisotropic damage model
enum dam_evolfunc {brittle = 1,  quasi_brittle= 2, quadbezier=3}; 
const enumstr dam_evolfuncstr[] = {{"brittle",1}, {"quasi_brittle",2}, {"quadbezier",3}}; 
const kwdset dam_evolfunc_kwdset(sizeof(dam_evolfuncstr)/sizeof(*dam_evolfuncstr), dam_evolfuncstr);



// aliases for switching on/off of fatigue in damage models
enum fatigue_flag {fatigue_off = 0, fatigue_on = 1};
const enumstr fatigue_flagstr[] = {{"fatigue_off",0}, {"fatigue_on",1}}; 
const kwdset fatigue_flag_kwdset(sizeof(fatigue_flagstr)/sizeof(*fatigue_flagstr), fatigue_flagstr);



// aliases for type of averaging in nonlocal material models
enum wavrg {avggamma = 1, avgepsp = 2};                  
const enumstr wavrgstr[] = {{"avggamma",1}, {"avgepsp",2}};
const kwdset wavrg_kwdset(sizeof(wavrgstr)/sizeof(*wavrgstr), wavrgstr);



// aliases for types of averaging functions in nonlocal material models
enum avrgf {parab = 1, cubic = 2, exponential = 3};
const enumstr avrgfstr[] = {{"parab",1}, {"cubic",2}, {"exponential",3}};
const kwdset avrgf_kwdset(sizeof(avrgfstr)/sizeof(*avrgfstr), avrgfstr);



// aliases for types of functions for modulus evolution in elasttime model
enum ym_evolfunc {b3law = 1,doublepwrlaw = 2,userfunc = 3};
const enumstr ym_evolfuncstr[] = {{"b3law",1}, {"doublepwrlaw",2}, {"userfunc",3}};
const kwdset ym_evolfunc_kwdset(sizeof(ym_evolfuncstr)/sizeof(*ym_evolfuncstr), ym_evolfuncstr);



// aliases for types of shrinkage law 
enum tshrlaw {shr_measured = 1, beta_val = 2};
const enumstr tshrlawstr[] = {{"shr_measured",1}, {"beta_val",2}};
const kwdset tshrlaw_kwdset(sizeof(tshrlawstr)/sizeof(*tshrlawstr), tshrlawstr);



// aliases for format of graphics output
enum graphfmt {grfmt_no=0, grfmt_open_dx=1, grfmt_femcad=2, grfmt_gid=3, grfmt_gid_sep=4, grfmt_vtk=5,
               grfmt_gid_vtk=6, grfmt_gidsep_vtk=7, grfmt_sep_files_quant=8};
const enumstr graphfmtstr[] = {{"grfmt_no",0}, {"grfmt_open_dx",1}, {"grfmt_femcad",2}, 
                               {"grfmt_gid",3}, {"grfmt_gid_sep",4}, {"grfmt_vtk",5}, 
                               {"grfmt_gid_vtk",6}, {"grfmt_gidsep_vtk",7}, {"grfmt_sep_files_quant",8}};
const kwdset graphfmt_kwdset(sizeof(graphfmtstr)/sizeof(*graphfmtstr), graphfmtstr);



// aliases for type of printed unknown
enum prunk {pr_displ=1, pr_strains=2, pr_stresses=3, pr_forces=4, pr_react=5, 
            pr_stepid=6, pr_appload=7, pr_other=8, pr_time=9, pr_eigval=10, pr_macrostrain=11, pr_macrostress=12, pr_residual=13};
const enumstr prunkstr[] = {{"pr_displ",1}, {"pr_strains",2}, {"pr_stresses",3}, {"pr_forces",4}, 
                            {"pr_react",5}, {"pr_stepid",6}, {"pr_appload",7}, {"pr_other",8}, 
                            {"pr_time",9},  {"pr_eigval",10},{"pr_macrostrain",11},{"pr_macrostress",12}, {"pr_residual",13}};
const kwdset prunk_kwdset(sizeof(prunkstr)/sizeof(*prunkstr), prunkstr);



//aliases for time printing
enum timetypeprin {seconds=1,minutes=2,hours=3,days=4};
const enumstr timetypeprinstr[] = {{"seconds",1}, {"minutes",2}, {"hours",3}, {"days",4}};
const kwdset timetypeprin_kwdset(sizeof(timetypeprinstr)/sizeof(*timetypeprinstr), timetypeprinstr);



// aliases for types of timedependent load
enum dynload {timeindload = 1, seismicload=10, responsespectrum=11, timedepload=20};
const enumstr dynloadstr[] = {{"timeindload",1}, {"seismicload",10}, {"responsespectrum",11}, {"timedepload",20}};
const kwdset dynload_kwdset(sizeof(dynloadstr)/sizeof(*dynloadstr), dynloadstr);



// aliases for directions of dynamic load
enum dirdynload {xdir=1,ydir=2,zdir=3};
const enumstr dirdynloadstr[] = {{"xdir",1}, {"ydir",2}, {"zdir",3}};
const kwdset dirdynload_kwdset(sizeof(dirdynloadstr)/sizeof(*dirdynloadstr), dirdynloadstr);



// aliases types of eigenstrains
enum eigstraintype {eigstrain=1,tempstrain=2};
const enumstr eigstraintypestr[] = {{"eigstrain",1}, {"tempstrain",2}};
const kwdset eigstraintype_kwdset(sizeof(eigstraintypestr)/sizeof(*eigstraintypestr), eigstraintypestr);



// aliases for types of graph approximation in graphmat material models
enum graphtype { glinear = 0, gtable = 1, gfunc = 2, gfunc_ser = 3};
const enumstr graphtypestr[] = {{"glinear",0}, {"gtable",1}, {"gfunc",2}, {"gfunc_ser",3}};
const kwdset graphtype_kwdset(sizeof(graphtypestr)/sizeof(*graphtypestr), graphtypestr);



// aliases for different types of hardening function in plasticity material models
enum hardensoften {nohs=0,linearhs=1,sqrths=2,sinhs=3,sqrtsinhs=4,
                   anlinearhs=101,ansqrths=102,ansinhs=103,ansqrtsinhs=104,
                   plstrainnorm=201,limplstrainnorm=202};
const enumstr hardensoftenstr[] = {{"nohs",0}, {"linearhs",1}, {"sqrths",2}, {"sinhs",3}, {"sqrtsinhs",4},
                                   {"anlinearhs",101}, {"ansqrths",102}, {"ansinhs",103}, {"ansqrtsinhs",104},
                                   {"plstrainnorm",201}, {"limplstrainnorm",202}};
const kwdset hardensoften_kwdset(sizeof(hardensoftenstr)/sizeof(*hardensoftenstr), hardensoftenstr);



// aliases for input variable to hardening/softening
enum hsinputvar {plstrnorm=1,consistpar=2,limplstrnorm=11,limconsistpar=12,limanplstrnorm=21,limanconsistpar=22};
const enumstr hsinputvarstr[] = {{"plstrnorm",1}, {"consistpar",2}, {"limplstrnorm",11}, {"limconsistpar",12},
                                 {"limanplstrnorm",21}, {"limanconsistpar",22}};
const kwdset hsinputvar_kwdset(sizeof(hsinputvarstr)/sizeof(*hsinputvarstr), hsinputvarstr);



// aliases for ???!!! probably not used - candidate for removal
enum strestate {notdef=0,tension=1,compression=2};
const enumstr strestatestr[] = {{"notdef",0}, {"tension",1}, {"compression",2}};
const kwdset strestate_kwdset(sizeof(strestatestr)/sizeof(*strestatestr), strestatestr);



// aliases for type of pore pressure computation
enum pore_press_comp {noporep=0, const_all=1, nonmechq=2, var_press=3};
const enumstr pore_press_compstr[] = {{"noporep",0}, {"const_all",1}, {"nonmechq",2}, {"var_press",3}};
const kwdset pore_press_comp_kwdset(sizeof(pore_press_compstr)/sizeof(*pore_press_compstr), pore_press_compstr);



// aliases for type of coupled
enum pore_press_coup {part_coup_pp=1, fully_coup_pp=2};
const enumstr pore_press_coupstr[] = {{"part_coup_pp",1}, {"fully_coup_pp",2}};
const kwdset pore_press_coup_kwdset(sizeof(pore_press_coupstr)/sizeof(*pore_press_coupstr), pore_press_coupstr);



// aliases for integrated quantity (used in elem_integration)
enum integratedquant {locstress=1,nonlocstress=2,stressincr=3,eigstress=4,penergydens=5};
const enumstr integratedquantstr[] = {{"locstress",1}, {"nonlocstress",2}, {"stressincr",3}, {"eigstress",4}, {"penergydens",5}};
const kwdset integratedquant_kwdset(sizeof(integratedquantstr)/sizeof(*integratedquantstr), integratedquantstr);

/// aliases for definition of axis (used in gradual construction for initial displacements given by rotation about axis)
enum axisdeftype {noaxdef=0, ptx2=1, nod_pt=2, nodx2=3};
const enumstr axisdeftypestr[] = {{"noaxdef",0}, {"ptx2",1}, {"nod_pt",2}, {"nodx2",3}};
const kwdset axisdeftype_kwdset(sizeof(axisdeftypestr)/sizeof(*axisdeftypestr), axisdeftypestr);



// aliases for the type of residual vector norm computation in the Newton-Raphson method
enum resnormt {rel_load_norm=1, rel_react_norm=2, rel_loadreact_norm=3, absol_norm=4};
const enumstr resnormtstr[] = {{"rel_load_norm",1}, {"rel_react_norm",2}, {"rel_loadreact_norm",3}, {"absol_norm",4}};
const kwdset resnormt_kwdset(sizeof(resnormtstr)/sizeof(*resnormtstr), resnormtstr);



// aliases for type of mechanical quantities
enum mechquant {
  nomech_q=0,
  displ_q=1,         // displacement vector (r)
  strain_q=2,        // strain tensor (eps)
  stress_q=3,        // stress tensor (sig)
  other_q=4,         // other array values as a vector (other)
  react_q=5,         // reactions vector at nodes  (R)              
  load_vect_q=6,      // load vector at nodes (f_l)
  int_force_q=7,     // internal force vector (f_i)
  resid_vect_q=8,     // residual vector (f_r)
  tempr_strain_q=9,  // temperature strain tensor (epst)
  eig_strain_q=10,   // eigenstrain tensor (eps0)
  eig_stress_q=11,   // eigenstress tensor (sig0)
  time_q=12,         // time (constant scalar quantity at whole domain)
  step_id_q=13,      // step id (constant scalar quantity at whole domain)
  load_fact_q=14,    // load factor in nonlinear statics problem type (constant scalar quantity at whole domain)
  eigval_q=15,       // eigen values in eigenvalue problem type (constant scalar quantity at whole domain)
  macrostrain_q=16,  // macrostrain (constant tensor quantity at whole domain) (Sig)
  macrostress_q=17,  // macrostress (constant tensor quantity at whole domain)  (Eps)
  nonmech_q=18,      // nonmechanical quantities at ip
  first_inv=19,       // A11 + A22 + A33 (I_1)
  second_inv=20,      // A11*A22 + A11*A33 + A22*A33 - A12^2 - A13^2 - A23^2 (I_2)
  third_inv=21,       // det|A| (I_3)
  tensor_norm=22,     // ||A|| = sqrt(a_ij*a_ij) (||.||)
  tensdeviator=23,    // D_ij = A_ij - delta_ij*A_kk/3 (d)
  strain_vol=24,      // eps_v = eps_x + eps_y + eps_z (eps_v)
  mean_stress=25,     // sig_m = (sig_x + sig_y + sig_z)/3 (sig_m)
  j2inv=26,           // negative value of the second invariant of stress deviator, i.e. J2 = 1/2 s_ij s_ij (J2)
  von_mises_stress=27,// sig_eff = sqrt(3*J2)  (sig_vM)              
  strain_deviator=28, // e_ij = eps_ij - delta_ij*eps_v/3  (e)
  stress_deviator=29, // s_ij = sig_ij - delta_ij*sig_m    (s)
  strain_pl=30,      // plastic strain tensor              (epsp)
  cons_param=31,     // consistency parameter gamma
  damage_scal=32,    // scalar damage (omega)
  damaget_scal=33,   // scalar damage (omega_t) in tension
  damagec_scal=34,   // scalar damage (omega_c) in compression
  damage_tens=35,    // damage tensor (Omega)
  damaget_tens=36,   // damage tensor (Omega_t) in tension
  damagec_tens=37,   // damage tensor (Omega_c) in compression
  eps_x=38, eps_y=39, eps_z=40, gamma_yz=41, gamma_xz=42, gamma_xy=43, eps_yz=44, eps_xz=45, eps_xy=46, // strain components as scalars
  sig_x=47, sig_y=48, sig_z=49, tau_yz=50, tau_xz=51, tau_xy=52, // stress components as scalars
  epsp_x=53, epsp_y=54, epsp_z=55, gammap_yz=56, gammap_xz=57, gammap_xy=58, epsp_yz=59, epsp_xz=60, epsp_xy=61  // components of plastic strains considered as scalars
};



const enumstr mechquantstr[] = {
  {"nomech_q", 0},
  {"displ_q", 1},         // displacement vector
  {"strain_q", 2},        // strain tensor
  {"stress_q", 3},        // stress tensor
  {"other_q", 4},         // other array values as scalar
  {"react_q", 5},         // reactions vector at nodes                
  {"load_vect_q", 6},     // load vector at nodes
  {"int_force_q", 7},     // internal force vector
  {"resid_vect_q", 8},     // residual vector
  {"tempr_strain_q", 9},  // temperature strain tensor
  {"eig_strain_q", 10},   // eigenstrain tensor
  {"eig_stress_q", 11},   // eigenstress tensor
  {"time_q", 12},         // time (constant scalar quantity at whole domain)
  {"step_id_q", 13},      // step id (constant scalar quantity at whole domain)
  {"load_fact_q", 14},    // load factor in nonlinear statics problem type (constant scalar quantity at whole domain)
  {"eigval_q", 15},       // eigen values in eigenvalue problem type (constant scalar quantity at whole domain)
  {"macrostrain_q", 16},  // macrostrain (constant tensor quantity at whole domain)
  {"macrostress_q", 17},  // macrostress (constant tensor quantity at whole domain)
  {"nonmech_q", 18},      // nonmechanical quantities at ip as a vector
  {"first_inv", 19},       // A11 + A22 + A33
  {"second_inv", 20},      // A11*A22 + A11*A33 + A22*A33 - A12^2 - A13^2 - A23^2
  {"third_inv", 21},       // det|A|
  {"tensor_norm", 22},     // ||A|| = sqrt(a_ij*a_ij)
  {"tensdeviator", 23},    // D_ij = A_ij - delta_ij*A_kk/3 
  {"strain_vol", 24},      // eps_v = eps_x + eps_y + eps_z
  {"mean_stress", 25},     // sig_m = (sig_x + sig_y + sig_z)/3
  {"j2inv", 26},           // negative value of the second invariant of stress deviator}, i.e. J2 = 1/2 s_ij s_ij
  {"von_mises_stress", 27},// sig_eff = sqrt(3*J2)                
  {"strain_deviator", 28}, // e_ij = eps_ij - delta_ij*eps_v/3
  {"stress_deviator", 29}, // s_ij = sig_ij - delta_ij*sig_m
  {"strain_pl", 30},       // plastic strain tensor
  {"cons_param", 31},      // consistency parameter gamma
  {"damage_scal", 32},     // scalar damage (omega)
  {"damaget_scal", 33},    // scalar damage (omega_t) in tension
  {"damagec_scal", 34},    // scalar damage (omega_c) in compression
  {"damage_tens", 35},     // damage tensor (Omega)
  {"damaget_tens", 36},    // damage tensor (Omega_t) in tension
  {"damagec_tens", 37},    // damage tensor (Omega_c) in compression
  {"eps_x", 38},           // strain components as scalars
  {"eps_y", 39},
  {"eps_z", 40},
  {"gamma_yz", 41},
  {"gamma_xz", 42},
  {"gamma_xy", 43},
  {"eps_yz", 44},
  {"eps_xz", 45},
  {"eps_xy", 46},
  {"sig_x", 47},           // stress components as scalars
  {"sig_y", 48},
  {"sig_z", 49},
  {"tau_yz", 50},
  {"tau_xz", 51},
  {"tau_xy", 52},
  {"epsp_x", 53},          // components of plastic strains considered as scalars
  {"epsp_y", 54},
  {"epsp_z", 55},
  {"gammap_yz", 56},
  {"gammap_xz", 57},
  {"gammap_xy", 58},
  {"epsp_yz", 59},
  {"epsp_xz", 60},
  {"epsp_xy", 61}
};
const kwdset mechquant_kwdset(sizeof(mechquantstr)/sizeof(*mechquantstr), mechquantstr);
#endif
