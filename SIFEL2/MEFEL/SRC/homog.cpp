#include "homog.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "element.h"
#include "intpoints.h"
#include "globmat.h"
#include "gtopology.h"
#include "math.h"

/*
routines to average strain and stress over the domain of hexahedral elements ONLY!!!
other elements will not work
Mp->homog==1 effective or apparent elastic modulus
Mp->homog==2 simulation of autogenous shrinkage, fixed lid problem
smilauer@cml.fsv.cvut.cz, 16.1.2006
*/

//function to cycle index in a 3 dimensional array
int C012(int x){
  if (x<0){
    x+=3;
  }
  else if(x>2){
    x-=3;
  }
  return x;
}

void homogenization (FILE */*out*/, long lcid)
{

  if (Mp->homog!=3){

  //for homogenization , eigenstrains are present in the strains
  //stresses are already calculated without eigenstrains

  //isotropic 3D material
  //compute_req_val (lcid);//add eigenstrains;
  
  long i, j, k, ncomp, id, ipp, tnipe, ndofn, counter;
  double *stress_av, *strain_av;
  double *stress_max, *strain_max;
  double *r;
  double x_sum, y_sum, z_max, z_sum;
  double stress_eq_max=0;
  double mu[3], E[9], nu[9];
  double G_en, nu_en, E_en, W_en, help_a, help_b;
  FILE *propfile;
  char filename[]="Elas_1.out";

  ncomp=6;//assume 3D elasticity

  stress_av = new double [ncomp]; //contains average stresses
  strain_av = new double [ncomp]; //contains average strains
  
  stress_max = new double [ncomp]; //contains maximum stresses
  strain_max = new double [ncomp]; //contains maximum strains

  if((propfile=fopen(filename, "a"))==NULL){
    printf("error in opening %s file for writing\n", filename);
    perror("");
  }
  
  for(i=0;i<ncomp;i++){
    stress_av[i]=0;
    strain_av[i]=0;
    stress_max[i]=0;
    strain_max[i]=0;
  }


  //store strains and stresses from all integration points
  for (i=0; i<Mt->ne; i++){
    ipp = Mt->elements[i].ipp[0][0];
    tnipe = Mt->give_totnip(i);
    ncomp = Mm->ip[ipp].ncompstr;
    id = lcid*ncomp;
    for (j=0; j<tnipe; j++){
      for (k=0; k<ncomp; k++){
	//substract eigenstrains
	help_a=Mm->ip[ipp+j].stress[id+k];
	help_b=Mm->ip[ipp+j].strain[id+k]-Mm->eigstrains[j][k];

	stress_av[k]+=help_a;
	strain_av[k]+=help_b;

	//strain_av[k]-=Mm->eigstrains[j][k];
	//strain_av[k]+=Mm->ip[ipp+j].strain[id+k];
	//store maximum stresses and strains 
	if(help_a<stress_max[k]){
	  stress_max[k]=help_a;
	}
	if(help_b<strain_max[k]){
	  strain_max[k]=help_b;
	}
      }

      help_a=sqrt(2.)/2.*
	sqrt((Mm->ip[ipp+j].stress[id+0]-Mm->ip[ipp+j].stress[id+1])*(Mm->ip[ipp+j].stress[id+0]-Mm->ip[ipp+j].stress[id+1])+
	     (Mm->ip[ipp+j].stress[id+1]-Mm->ip[ipp+j].stress[id+2])*(Mm->ip[ipp+j].stress[id+1]-Mm->ip[ipp+j].stress[id+2])+
	     (Mm->ip[ipp+j].stress[id+2]-Mm->ip[ipp+j].stress[id+0])*(Mm->ip[ipp+j].stress[id+2]-Mm->ip[ipp+j].stress[id+0])+
	     6*Mm->ip[ipp+j].stress[id+3]*Mm->ip[ipp+j].stress[id+3]+
	     6*Mm->ip[ipp+j].stress[id+4]*Mm->ip[ipp+j].stress[id+4]+
	     6*Mm->ip[ipp+j].stress[id+5]*Mm->ip[ipp+j].stress[id+5]);
      
      if(help_a>stress_eq_max){
	stress_eq_max=help_a;
      }

    }
  }
  
  //AVERAGE ENGINEERING STRAIN & STRESS
  for(i=0;i<ncomp;i++){
    //perform averaging
    strain_av[i]/=(-Mt->ne*tnipe);
    stress_av[i]/=(-Mt->ne*tnipe);
  }
  
  //calculate E and nu from the energetic equivalence - only 3D isotropic material
  if (Mp->homog==1){//case of averaging 

 
    G_en=(strain_av[3]*stress_av[3]+strain_av[4]*stress_av[4]+
	  strain_av[5]*stress_av[5])/
      (strain_av[3]*strain_av[3]+strain_av[4]*strain_av[4]+
       strain_av[5]*strain_av[5]);
  
    help_a=(strain_av[0]*strain_av[0]+strain_av[1]*strain_av[1]+
	    strain_av[2]*strain_av[2]);
  
    help_b=2*(strain_av[0]*strain_av[1]+strain_av[0]*strain_av[2]+
	      strain_av[1]*strain_av[2]);
  
  
    W_en=0.5*(strain_av[0]*stress_av[0]+strain_av[1]*stress_av[1]+
	      strain_av[2]*stress_av[2]);
  
    nu_en=(W_en-G_en*help_a)/(G_en*(help_b-help_a)+2*W_en);
    E_en=2*G_en*(1+nu_en);
  
    //print elastic constants E, nu, G=mu, k
    fprintf(propfile, "ENER. %f %f %f %f %f ", E_en, nu_en, G_en, E_en/(3.*(1.-2*nu_en)), W_en);
  
    fprintf(propfile, "STRE ");
    for(i=0;i<ncomp;i++){
      fprintf(propfile,"%f ", stress_av[i]);
    }
  
    fprintf(propfile,"STRN ");
    for(i=0;i<ncomp;i++){
      fprintf(propfile,"%f ", strain_av[i]);
    }

    //equivalent M-H-H stress 
    fprintf(propfile,"MHH ");
    fprintf(propfile,"%f ", stress_eq_max);


    fprintf(propfile,"MAX_STRE ");
    for(i=0;i<ncomp;i++){
      fprintf(propfile,"%f ", stress_max[i]);
    }

    fprintf(propfile,"MAX_STRN ");
    for(i=0;i<ncomp;i++){
      fprintf(propfile,"%f ", strain_max[i]);
    }
  
    //COMPUTE ELASTIC SHEAR MODULI in x,y,z direction
    mu[0]=stress_av[3]/strain_av[3];
    mu[1]=stress_av[4]/strain_av[4];
    mu[2]=stress_av[5]/strain_av[5];
  
    //COMPUTE ELASTIC YOUNG'S MODULI in x,y,z direction
    fprintf(propfile, "E ");
    for (i=0; i<=2; i++){//x,y,z
      for (j=0; j<=2; j++){//components in one direction
	E[3*i+j]=(stress_av[C012(i+0)]+stress_av[C012(i+1)]+stress_av[C012(i+2)])/
	  (strain_av[C012(i+0)]+(stress_av[C012(i+1)]+stress_av[C012(i+2)])/(2*mu[j]));
	fprintf(propfile, "%f ", E[3*i+j]);
      }
    }
  
    //COMPUTE POISSON's AND SHEAR MODULI in x,y,z direction
    fprintf(propfile, "Nu ");
    for (i=0; i<=2; i++){//x,y,z
      for (j=0; j<=2; j++){//components in one direction
	nu[3*i+j]=(E[3*i+j]-2*mu[i])/(2*mu[i]);
	fprintf(propfile, "%f ", nu[3*i+j]); 
      }
    }
  
    fprintf(propfile, "\n");

  }

  /*
  //average top displacements, values are too much overestimated when all phase taken into account, only solid voxels are considered 
  else if(Mp->homog==2){
    double z_av;

    //find maximum z coordinate of all nodes
    z_max=0.;
    for (i=0; i<Mt->nn; i++){
      if(Gtm->gnodes[i].z>z_max){
	z_max=Gtm->gnodes[i].z;
      }
    }
    
    z_sum=0.;
    counter=0;
    for (i=0; i<Mt->nn; i++){
      ndofn = Mt->give_ndofn(i);
      r = new double[ndofn];
      noddispl(0, r, i);//displacements of load case 0
      if(z_max==Gtm->gnodes[i].z){
	//printf("node %ld at %f %f %f r_z %e\n", i, Gtm->gnodes[i].x, Gtm->gnodes[i].y, Gtm->gnodes[i].z, r[2]);
	z_sum+=r[2];
	counter++;
	delete [] r; 
      }
    } 
    help_a=counter;
    z_av=z_sum/counter;
    
    //now consider only displacements that are smaller than avegage (aim = exclude porosity)
    counter=0;
    for (i=0; i<Mt->nn; i++){
      ndofn = Mt->give_ndofn(i);
      r = new double[ndofn];
      noddispl(0, r, i);//displacements of load case 0
	if((z_max==Gtm->gnodes[i].z) && (r[2]<z_av)){
	  z_sum+=r[2];
	  counter++;
	}
	delete [] r;
    }
    
    printf("\ntop_nodes %.0f / %d z_av %f z_av_sel %f\n", help_a, counter, z_av, z_sum/counter);
      
    //printf("top_nodes %ld z_sum %e average_eps %e\n", counter, z_sum, z_sum/counter/z_max);
    //fprintf(propfile, "top_nodes %ld z_sum %e average_eps %e\n", counter, z_sum, z_sum/counter/z_max);
  }
*/  

  
  else if(Mp->homog==2){//report stress and strain components from various phases - autogeneous shrinkage
    
    //find maximum z coordinate of all nodes - z is bigger due to the LID!!!!
    z_max=0.;
    x_sum=0.;
    y_sum=0.;
    z_sum=0.;
    for (i=0; i<Mt->nn; i++){
      if(Gtm->gnodes[i].z>z_max){
	z_max=Gtm->gnodes[i].z;
      }
    }
    
    printf("\n ");
    
    for(i=0; i<Mt->nn; i++){
      if(z_max==Gtm->gnodes[i].z){
	ndofn = Mt->give_ndofn(i);
	r = new double[ndofn];
	noddispl (lcid, r, i);
	x_sum+=r[0];
	y_sum+=r[1];
	z_sum+=r[2];
	//printf("%f\n", Gtm->gnodes[i].x,Gtm->gnodes[i].y,Gtm->gnodes[i].z,r[2]);
	delete [] r;
      }
    }
 
    x_sum/=(z_max*z_max);
    y_sum/=(z_max*z_max);
    z_sum/=(z_max*z_max);
   
    fprintf(propfile, "z_strn %f x_strn %f y_strn %f av_displ_z %f ", z_sum/(z_max-1), x_sum/(z_max-1), y_sum/(z_max-1), z_sum); 

    for(i=0;i<ncomp;i++){
      stress_av[i]=0;
      strain_av[i]=0;
      stress_max[i]=0;
      strain_max[i]=0;
    }

    //***************POROSITY*******************
    //calculate stress and strain for voxel, not what happens in water (equilibrium), i.e.
    //tension in water = compression in voxel
    for (k=0; k<ncomp; k++){
      stress_av[k]=0;
      strain_av[k]=0;
    }
    counter=0;
   
    for (i=0; i<Mt->ne; i++){
      //empty porosity (26-1), water-filled porosity (phase 28-1), Vit Smilauer, 10.1.2006
      if (Mt->elements[i].idm[0]==(28-1)){
	ipp = Mt->elements[i].ipp[0][0];
	tnipe = Mt->give_totnip(i);
	ncomp = Mm->ip[ipp].ncompstr;
	id = lcid*ncomp;
	for (j=0; j<tnipe; j++){
	  counter++;
	  for (k=0; k<ncomp; k++){
	    //do not substract eigenstrains, use values with eigenstrains
	    help_a=Mm->ip[ipp+j].stress[id+k];
	    help_b=Mm->ip[ipp+j].strain[id+k];//-Mm->eigstrains[j][k];
	    
	    stress_av[k]+=help_a;
	    strain_av[k]+=help_b;
	  }
	}
      }
    } 
  
    //AVERAGE ENGINEERING STRAIN & STRESS
    for(i=0;i<ncomp;i++){
      //perform averaging
      strain_av[i]/=(counter);
      stress_av[i]/=(counter);
    }
  
    fprintf(propfile,"SIG_POR "); 

    //average stress from selected integration points
    for(i=0;i<ncomp;i++){
      fprintf(propfile,"%f ", stress_av[i]);
    }

    fprintf(propfile,"EPS_POR ");     
    //average strain from selected integration points
    for(i=0;i<ncomp;i++){
      fprintf(propfile,"%f ", strain_av[i]);
    }

   /************only SOLIDS, themselves or in RVE****************/
    for (k=0; k<ncomp; k++){
      stress_av[k]=0;
      strain_av[k]=0;
    }
    counter=0;
    
    for (i=0; i<Mt->ne; i++){
      //empty porosity (26-1), water-filled porosity (phase 28-1), stiff lid (phase(29-1)
      if(Mt->elements[i].idm[0]!=(28-1) && Mt->elements[i].idm[0]!=(26-1) && Mt->elements[i].idm[0]!=(29-1)){
	ipp = Mt->elements[i].ipp[0][0];
	tnipe = Mt->give_totnip(i);
	ncomp = Mm->ip[ipp].ncompstr;
	id = lcid*ncomp;
	for (j=0; j<tnipe; j++){
	  counter++;
	  for (k=0; k<ncomp; k++){
	    //do not substract eigenstrains - makes no difference
	    help_a=Mm->ip[ipp+j].stress[id+k];
	    help_b=Mm->ip[ipp+j].strain[id+k];//-Mm->eigstrains[j][k];
	    stress_av[k]+=help_a;
	    strain_av[k]+=help_b;
	  }
	}
      }
    } 
  
    //averaging in RVE
    fprintf(propfile,"SIG_SOL_RVE ");
    //average stress from selected integration points - 8 integration points on 1 element
    for(i=0;i<ncomp;i++){
      fprintf(propfile,"%f ", stress_av[i]/(z_max-1.)/(z_max-1.)/(z_max-1.)/8.);
    }
    
    fprintf(propfile,"EPS_SOL_RVE ");     
    //average strain from selected integration points
    for(i=0;i<ncomp;i++){
      fprintf(propfile,"%f ", strain_av[i]/(z_max-1.)/(z_max-1.)/(z_max-1.)/8.);
    }

    fprintf(propfile,"SIG_SOL "); 
    //average stress from selected integration points
    for(i=0;i<ncomp;i++){
      fprintf(propfile,"%f ", stress_av[i]/counter);
    }

    fprintf(propfile,"EPS_SOL ");     
    //average strain from selected integration points
    for(i=0;i<ncomp;i++){
      fprintf(propfile,"%f ", strain_av[i]/counter);
    }
 

    fprintf(propfile,"\n"); 

  }//end of elseif

  fclose(propfile);

  delete [] stress_av;
  delete [] strain_av;
  delete [] stress_max;
  delete [] strain_max;

  }//  end of the statement if (Mp->homog!=3)
}
