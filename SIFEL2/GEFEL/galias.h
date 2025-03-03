#ifndef GALIAS_H
#define GALIAS_H
#include "kwdset.h"


// aliases for element types
enum gtypel {noel=0, isolinear1d=1, isoquadratic1d=2, trianglelinear=3, trianglequadratic=4, isolinear2d=5, isoquadratic2d=6, 
             tetrahedronlinear=7, tetrahedronquadratic=8, pyramidelinear=9, pyramidequadratic=10, wedgelinear=11,
             wedgequadratic=12, isolinear3d=13, isoquadratic3d=14, isocubic2d=15, all_element_types=16};
const enumstr gtypelstr[] = {{"noel",0}, {"isolinear1d", 1}, {"isoquadratic1d", 2}, {"trianglelinear", 3}, {"trianglequadratic", 4}, 
                             {"isolinear2d", 5}, {"isoquadratic2d", 6}, {"tetrahedronlinear", 7}, {"tetrahedronquadratic", 8}, 
                             {"pyramidelinear", 9}, {"pyramidequadratic", 10}, {"wedgelinear", 11}, {"wedgequadratic", 12},
                             {"isolinear3d", 13}, {"isoquadratic3d", 14}, {"isocubic2d", 15}};
const kwdset gtypel_kwdset(sizeof(gtypelstr)/sizeof(*gtypelstr), gtypelstr);



// aliases for basis functions
enum gtypebf {nobf=0, bar_linear01=1, bar_linear11=2, bar_quadratic01=3, bar_quadratic11=4,
	      beam2D_cubic01=11,beam2D_cubic11=12,
	      quadrilateral_bilinear=21,quadrilateral_biquadratic=22,quadrilateral_rotdofs=23,
	      triangle_linear=31,triangle_quadratic=32,triangle_rotdof=33,
	      cct=41,batoz_plate=42,
	      hexahedral_trilinear=51,hexahedral_triquadratic=52,
	      tetrahedral_linear=61,tetrahedral_quadratic=62};
const enumstr gtypebfstr[] = {{"nobf",0}, {"bar_linear01", 1}, {"bar_linear11", 2}, {"bar_quadratic01", 3}, {"barquadratic", 4}, 
			      {"beam2d_cubic01", 11}, {"beam2d_cubic11", 12},
			      {"quadrilateral_bilinear", 21}, {"quadrilateral_biquadratic", 22}, {"quadrilateral_rotdof", 23},
			      {"triangle_linear", 31}, {"triangle_quadratic", 32}, {"triangle_rotdof", 33},
			      {"cct", 41}, {"batoz_plate",42},
			      {"hexahedral_trilinear", 51}, {"hexahedral_triquadratic", 52},
			      {"tetrahedral_linear", 61}, {"tetrahedral_quadratic", 62}};
const kwdset gtypebf_kwdset(sizeof(gtypebfstr)/sizeof(*gtypebfstr), gtypebfstr);


// aliases for entity types in T3d
enum entityp {evertex=1, ecurve=2, esurface=3, eregion=4, epatch=5, eshell=6};
const enumstr entitypstr[] = {{"evertex", 1}, {"ecurve", 2}, {"esurface", 3}, {"eregion", 4}, {"epatch", 5}, {"eshell", 6}};
const kwdset entityp_kwdset(sizeof(entitypstr)/sizeof(*entitypstr), entitypstr);



// aliases for general enities
enum gentity {gvertex=1, gcurve=2, gsurface=3, gregion=4};
const enumstr gentitystr[] = {{"gvertex",1}, {"gcurve",2}, {"gsurface",3}, {"gregion",4}};
const kwdset gentity_kwdset(sizeof(gentitystr)/sizeof(*gentitystr), gentitystr);



//  aliases for mesh description
enum meshdescription {all_nodes=1,bound_nodes=2,neg_bound_nodes=3,glob_glued=4,metis=11,metis_elem=12};
const enumstr meshdescriptionstr[] = {{"all_nodes",1}, {"bound_nodes",2}, {"neg_bound_nodes",3}, {"glob_glued",4},
				      {"metis",11}, {"metis_elem",12}};
const kwdset meshdescription_kwdset(sizeof(meshdescriptionstr)/sizeof(*meshdescriptionstr), meshdescriptionstr);



// aliases for matrix storage in memory
enum storagetype {without_matrix=0,
                  dense_matrix=1,skyline_matrix=2,double_skyline=3,
                  compressed_rows=10,symm_comp_rows=11,compressed_columns=12,symm_comp_columns=13,element_matrices=40,diag_mat=50,
                  lapack_stor=100,petsc=120,spdirect_stor_scr=140,spdirect_stor_cr=141};
const enumstr storagetypestr[] = {{"without_matrix",0}, {"dense_matrix",1}, {"skyline_matrix", 2},
                                  {"double_skyline",3}, {"compressed_rows",10},{"symm_comp_rows", 11},
                                  {"compressed_columns",12},{"symm_comp_columns", 13},
                                  {"element_matrices", 40}, {"diag_mat",50}, 
				  {"lapack_stor", 100}, {"petsc", 120}, 
                                  {"spdirect_stor_scr",140}, {"spdirect_stor_cr", 141}};
const kwdset storagetype_kwdset(sizeof(storagetypestr)/sizeof(*storagetypestr), storagetypestr);



// aliases for solvers of linear equation system
enum linsolvertype {no_solver=0,gauss_elim=1,ldl=2,lu=3,ll=4,ill=5,conden=11,cg=20,
		    jacobi=24, gauss_seidel=25,bicg=30,diagsolv=6,
                    lapack_sol=100,pardisolver=120,spdirldl=140,spdirlu=141,spdirll=142,
                    sschur=151,sfeti=161,sfetidef=171, saddle_point=181, permon_solver=201, cholmod_solver=240};
const enumstr linsolvertypestr[] = {{"no_solver",0}, {"gauss_elim",1}, {"ldl",2}, {"lu",3}, {"ll",4}, {"ill",5}, {"conden",11},
                                    {"cg",20}, {"jacobi",24}, {"gauss_seidel",25},
				    {"bicg",30}, {"lapack_sol",100}, {"pardisolver",120}, {"spdirldl",140},
                                    {"spdirlu",141}, {"spdirll",142}, {"diagsolv",6},
				    {"sschur",151}, {"sfeti",161}, {"sfetidef",171},{"saddle_point",181},{"permon_solver",201}, 
                                    {"cholmod_solver", 240}};
const kwdset linsolvertype_kwdset(sizeof(linsolvertypestr)/sizeof(*linsolvertypestr), linsolvertypestr);



//  aliases for type of FETI implementation
enum fetiimplem {no_impl=0,boolean_matrices=1,nonredundant=2,redundant=3};
const enumstr fetiimplemstr[] = {{"no_impl",0}, {"boolean_matrices",1}, {"nonredundant",2}, {"redundant",3}};
const kwdset fetiimplem_kwdset(sizeof(fetiimplemstr)/sizeof(*fetiimplemstr), fetiimplemstr);



//  aliases for type of solver of coarse problem
enum redsystsolver {master_sol=1,paral_sol=2};
const enumstr redsystsolverstr[] = {{"master_sol",1}, {"paral_sol",2}};
const kwdset redsystsolver_kwdset(sizeof(redsystsolverstr)/sizeof(*redsystsolverstr), redsystsolverstr);



// aliases for preconditioners for conjugate gradients
enum precondtype {noprecond=0,diagprec=1,ssorprec=5,incomdec=10,sparseindec=20,boss=101,feti_lumped=151,feti_dirichlet=152};
const enumstr precondtypestr[] = {{"noprecond",0}, {"diagprec",1}, {"ssorprec",5}, {"incomdec",10}, {"sparseindec",20}, {"boss",101},{"feti_lumped",151},{"feti_dirichlet",152}};
const kwdset precondtype_kwdset(sizeof(precondtypestr)/sizeof(*precondtypestr), precondtypestr);



// aliases for node renumbering
enum noderenumb {no_renumbering=0,cuthill_mckee=1,rev_cuthill_mckee=2,sloan=3};
const enumstr noderenumbstr[] = {{"no_renumbering",0}, {"cuthill_mckee",1}, {"rev_cuthill_mckee",2}, {"sloan",3}};
const kwdset noderenumb_kwdset(sizeof(noderenumbstr)/sizeof(*noderenumbstr), noderenumbstr);



// aliases for eigenvalues solver
enum eigensolver {no_eigsolver=0,inv_iteration=1,subspace_it_jacobi=4,subspace_it_gsortho=5,shifted_subspace_it_gsortho=6};
const enumstr eigensolverstr[] = {{"no_eigsolver",0}, {"inv_iteration",1}, {"subspace_it_jacobi",4}, {"subspace_it_gsortho",5}, {"shifted_subspace_it_gsortho",6}};
const kwdset eigensolver_kwdset(sizeof(eigensolverstr)/sizeof(*eigensolverstr), eigensolverstr);



// aliases for different types of general function definition
enum generalfunct {constant=0,pars=1,tab=2,pars_set=3,tab_file=4,tab_file2=5,multpv=10,itab=20,gfunc_set=30};
const enumstr generalfunctstr[] = {{"stat",0}, {"pars",1}, {"tab",2}, {"pars_set",3}, {"tab_file",4}, {"tab_file2",5}, {"multpv", 10}, {"itab",20}, {"gfunc_set",30}};
const kwdset generalfunct_kwdset(sizeof(generalfunctstr)/sizeof(*generalfunctstr), generalfunctstr);



// aliases for different types of selection
enum seltype {sel_no=0, sel_all=1, sel_range=2, sel_list=3, sel_period=4, 
              sel_realrange=5, sel_reallist=6, sel_mtx=7, sel_range_mtx=8, sel_range_vec=9,
              sel_realperiod=10, sel_impvalues=11, sel_vec=12, sel_prop=13, sel_impvallst=14, sel_sarray=15, sel_single=16};
const enumstr seltypestr[] = {{"sel_no",0}, {"sel_all",1}, {"sel_range",2}, {"sel_list",3}, 
                              {"sel_period",4}, {"sel_realrange",5}, {"sel_reallist",6}, 
			      {"sel_mtx",7}, {"sel_rangem",8}, {"sel_rangev",9}, 
                              {"sel_realperiod",10}, {"sel_impvalues",11}, {"sel_vec",12},
                              {"sel_prop",13}, {"sel_impvallst",14}, {"sel_sarray",15}, {"sel_single", 16}};
const kwdset seltype_kwdset(sizeof(seltypestr)/sizeof(*seltypestr), seltypestr);



// aliases for point specification
enum nodip {no_point=0, atnode=1, atip=2, atxyz=3};
const enumstr nodipstr[] = {{"no_point",0}, {"atnode",1}, {"atip",2}, {"atxyz",3}};
const kwdset nodip_kwdset(sizeof(nodipstr)/sizeof(*nodipstr), nodipstr);

// aliases for simple flag
enum flagsw {off=0, on=1};
const enumstr flagswstr[] = {{"off",0}, {"on",1}};
const kwdset flagsw_kwdset(sizeof(flagswstr)/sizeof(*flagswstr), flagswstr);



// aliases for simple answer type
enum answertype {no=0, yes=1};
const enumstr answertypestr[] = {{"no",0}, {"yes",1}};
const kwdset answertype_kwdset(sizeof(answertypestr)/sizeof(*answertypestr), answertypestr);



// aliases for header/footer flag
enum hflagsw {hf_off=0, header=1, footer=2};
const enumstr hflagswstr[] = {{"hf_off",0}, {"header",1}, {"footer",2}};
const kwdset hflagsw_kwdset(sizeof(hflagswstr)/sizeof(*hflagswstr), hflagswstr);

// aliases for general types of finite elements
enum gelemtype {noelem=0,
		linbar=1, quadbar=2,
                lintriag=21,quadtriag=22,linquad=25,quadquad=26,cubicquad=27,
                lintetra=41,quadtetra=42,linhexa=45,quadhexa=46,linwed=47,gen2d=50};
const enumstr gelemtypestr[] = {{"noelem",0},
				{"linbar",1},{"quadbar",2},
				{"lintriag",21},{"quadtriag",22},{"linquad",25},{"quadquad",26},{"cubicquad",27},
				{"lintetra",41},{"quadtetra",42},{"linhexa",45},{"quadhexa",46},{"linwed",47}, 
                                {"gen2d",50}};
const kwdset gelemtype_kwdset(sizeof(gelemtypestr)/sizeof(*gelemtypestr),gelemtypestr);


// aliases for object type (used in the class sel)
enum objtype {gnod=1, gelem=2};

//  aliases for mesh format
enum meshform {sifel=0, t3d=1};
const enumstr meshformstr[] = {{"sifel",0}, {"t3d",1}};
const kwdset meshform_kwdset(sizeof(meshformstr)/sizeof(*meshformstr), meshformstr);



//  aliases for type of time controller
enum timecontrtype {notct=-1, fixed=0, adaptive=1, adaptivemin=2, adaptivemax=3, adaptive_minmax=4};
const enumstr timecontrtypestr[] = {{"notct",-1}, {"fixed",0}, {"adaptive",1}, {"adaptivemin",2}, {"adaptivemax",3}, {"adaptive_minmax",4}};
const kwdset timecontrtype_kwdset(sizeof(timecontrtypestr)/sizeof(*timecontrtypestr), timecontrtypestr);


// aliases for harddisk backup type
enum hdbackuptype {nohdb=0, hdbs_single=1, hdbs_multiple=2, hdbr_single=3, hdbr_multiple=4, hdbrs_single=5, 
                   hdbrs_multiple=6, hdbr_nonloc=7, hdbs_nonloc=8, hdbrs_nonloc=9};
const enumstr hdbackuptypestr[] = {{"nohdb",0}, {"hdbs_single",1}, {"hdbs_multiple",2}, {"hdbr_single",3}, {"hdbr_multiple",4}, 
                                   {"hdbrs_single",5}, {"hdbrs_multiple",6}, {"hdbr_nonloc",7}, {"hdbs_nonloc",8}, 
                                   {"hdbrs_nonloc",9}};
const kwdset hdbackuptype_kwdset(sizeof(hdbackuptypestr)/sizeof(*hdbackuptypestr), hdbackuptypestr);


// aliases for harddisk backup file format
enum hdbackupfmttype {text=1, binary=2};
const enumstr hdbackupfmttypestr[] = {{"text",1}, {"binary",2}};
const kwdset hdbackupfmttype_kwdset(sizeof(hdbackupfmttypestr)/sizeof(*hdbackupfmttypestr), hdbackupfmttypestr);


// aliases for keyword processing switch
enum pkwd_sw {nokwd=0, kwd_pd=1, kwd_outd=2, kwd_pd_outd=3, kwd_full=4};


// aliases for edge types used in material interfaces, jumps, hemivariational inequalities, etc.
enum edgetype {noedge=0,jumps=1, mater=2};
const enumstr edgetypestr[] = {{"noedge",0}, {"jumps",1}, {"mater",2}};
const kwdset edgetype_kwdset(sizeof(edgetypestr)/sizeof(*edgetypestr), edgetypestr);

//  aliases for type of macro-micro problem correspondence
enum macromicrotype {elem_domain=1,aggregate_domain=2,material_aggregate_domain=3,eff_aggregate_domain=4};
const enumstr macromicrotypestr[] = {{"elem_domain",1}, {"aggregate_domain",2}, {"eff_aggregate_domain",3}, {"eff_aggregate_domain",4}};
const kwdset macromicrotype_kwdset(sizeof(macromicrotypestr)/sizeof(*macromicrotypestr),macromicrotypestr);



// aliases for non-mechanical quantity (MUST be numbered from 1 and one by one) 
enum nonmechquant {temperature=1, rel_hum=2, initial_temperature=3, water_pressure=4,
                   cap_pressure=5, saturation_degree=6, suction=7, pore_pressure=8, vol_moist_cont=9,
                   free_salt_concentr=10, gas_pore_pressure=11,volume_change=12,eff_pore_pressure=13};
const enumstr nonmechquantstr[] = {{"temperature",1}, {"rel_hum",2}, {"initial_temperature",3},
                                   {"water_pressure",4}, {"cap_pressure",5}, {"saturation_degree",6}, 
                                   {"suction",7}, {"pore_pressure",8}, {"vol_moist_cont",9},
                                   {"free_salt_concentr",10}, {"gas_pore_pressure", 11}, {"volume_change",12}, {"eff_pore_pressure",13}};
const kwdset nonmechquant_kwdset(sizeof(nonmechquantstr)/sizeof(*nonmechquantstr), nonmechquantstr);



// aliases for non-transport quantity (MUST be numbered from 1 and one by one) 
enum nontransquant {precons_press=1, mean_stress_eff=2, virgin_porosity=3, init_porosity=4,
                    scal_iso_damage=5,proc_zone_length=6,crack_width=7,saturation_deg=8,der_saturation_deg=9,
                    porosity=10,void_ratio=11,advect_vel_x=12,advect_vel_y=13,advect_vel_z=14,der_saturation_deg_dtemp=15,
                    strain_vol_rate=16,der_saturation_deg_depsv=17, bulk_modulus=18, mmean_stress=19};
const enumstr nontransquantstr[] = {{"precons_press",1}, {"mean_stress_eff",2}, {"virgin_porosity",3}, 
                                    {"init_porosity",4}, {"scal_iso_damage",5}, {"proc_zone_length",6}, 
				    {"crack_width",7}, {"saturation_deg",8}, {"der_saturation_deg",9},
                                    {"porosity",10}, {"void_ratio",11}, {"advect_vel_x",12}, {"advect_vel_y",13},
				    {"advect_vel_z",14}, {"der_saturation_deg_dtemp",15}, {"strain_vol_rate",16}, 
				    {"der_saturation_deg_depsv",17}, {"bulk_modulus",18}, {"mmean_stress",19}};
const kwdset nontransquant_kwdset(sizeof(nontransquantstr)/sizeof(*nontransquantstr), nontransquantstr);



// aliases for quantity mathematical representation
enum quantrep {undefq=0, scalq=1, vectq=2, tensq=3};
const enumstr quantrepstr[] = {{"undefq", 0}, {"scalq", 1}, {"vectq", 2}, {"tensq", 3}};
const kwdset quantrep_kwdset(sizeof(quantrepstr)/sizeof(*quantrepstr), quantrepstr);


// aliases for the second order tensor storage notation
enum tensqnot {undeftn=0, voigtred=1, voigtfull=2, fullmtx=3};
const enumstr tensqnotstr[] {{"undeftn",0}, {"voigtred",1}, {"voigtfull",2}, {"fullmtx",3}};
const kwdset tensqnot_kwdset(sizeof(tensqnotstr)/sizeof(*tensqnotstr), tensqnotstr);


// aliases for the definition of local coordinate systems 
enum tmat_type {no_tr=0, gen_tr=1, cyl_tr=2, gfgen_tr=3}; 
const enumstr tmat_typestr[] {{"no_tr",0}, {"gen_tr",1}, {"cyl_tr",2}, {"gfgen_tr",3}}; 
const kwdset tmat_type_kwdset(sizeof(tmat_typestr)/sizeof(*tmat_typestr), tmat_typestr);

// aliases for format of graphics output
enum resfilefmt {resfmt_no=0, resfmt_open_dx=1, resfmt_femcad=2, resfmt_gid=3, resfmt_vtk=4,
                 resfmt_plain_out=5, resfmt_diag_dat=6};
const enumstr resfilefmtstr[] = {{"resfmt_no",0}, {"resfmt_open_dx",1}, {"resfmt_femcad",2},
                                 {"resfmt_gid",3}, {"resfmt_vtk",4}, {"resfmt_plain_out",5},
                                 {"resfmt_diag_dat",6}};
const kwdset resfilefmt_kwdset(sizeof(resfilefmtstr)/sizeof(*resfilefmtstr), resfilefmtstr);

enum timespectype {no_tspec=0, tspec_per, tspec_udt, tspec_tabper, tspec_selid};
const enumstr timespectypestr[] = {{"no_tspec",0}/*, tspec_per, tspec_udt, tspec_tabper, tspec_selid*/};
const kwdset timespectype_kwdset(sizeof(timespectypestr)/sizeof(*timespectypestr), timespectypestr);
#endif
