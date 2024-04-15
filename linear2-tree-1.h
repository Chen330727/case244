#include "heights.h"
#include "fractions.h"

#include "./linear2-tree-2.h"
#include "./vof2front-advanced-copy.h"

extern face vector modphase0,modphase1;
extern vector hhh;
extern int globali,outstep,level_interface;
extern scalar css_test3_n;
extern scalar deltac;
extern scalar T_solid;

extern scalar css_test2;
extern face vector fss_test2;

extern scalar css_test;
extern face vector fss_test;

//extern scalar topo_mask_g;
#define is_phase(v) (v>=0.5)
#define EPS 0.0000000001
 #define  is_boundary_box2(cell) is_boundary(cell)

foreach_dimension()
static inline void face_rejection2_x (Point point, scalar s)
{
  vector f = s.v;
  for (int i = 0; i <= 1; i++) {
    f.x[i] = 1.0;
#if dimension >= 2  
    f.x[i,1] = 1.0;
#endif
#if dimension >= 3
    for (int j = 0; j <= 1; j++)
      f.x[i,j,1] = 1.0;
#endif
  }
}

#include "myc2d2.h"
 extern  vector smallmodl,bigmodl,smallmodg,bigmodg;


struct Threephases {
  scalar ff4;
  face vector modphase1;
  face vector modphase0;
};

// void get_modphase01_2(struct Threephases q){
  void get_modphase01_2(){


//check mod//initial
foreach(){
    foreach_dimension(){
      smallmodl.x[]=-10.0;
      bigmodl.x[]=-10.0;
      smallmodg.x[]=-10.0;
      bigmodg.x[]=-10.0;

    }
}
///////////////

    vector nn4[];
    scalar alpha4[];
    scalar ff4[];//=q.ff4;
    // face vector modphase1=q.modphase1;
    // face vector modphase0=q.modphase0;
    // vector hhh[];
    foreach(){
       ff4[]=ff[];
    }   
    heights(ff4,hhh);
    boundary({hhh});
    ff4.restriction = restriction_volume_average;
    ff4.refine = ff4.prolongation = fraction_refine;
    restriction({ff4});

//taken from restruction
foreach() {

    if (ff4[] <= 0. || ff4[] >= 1.) {
      alpha4[] = 0.;
      foreach_dimension()
	       nn4.x[] = 0.;
      }else{
      coord m = interface_normal3 (point, ff4, hhh);
      foreach_dimension()
	        nn4.x[] = m.x;
      alpha4[] = plane_alpha (ff4[], m);
    }
  }

#if TREE
  foreach_dimension()
    nn4.x.refine = nn4.x.prolongation = refine_injection;

  alpha4.n = nn4;
  alpha4.refine = alpha4.prolongation = alpha_refine;
#endif
  

  //double lim_cut = 1e-12;
  double lim_cut = 1e-2;
  restriction({ff4,nn4,alpha4});
  face vector complete1[];
  // for(int l=0;l<=depth();l++){
  //      foreach_level(l){
  //         foreach_dimension(){
  //             complete1.x[] = 0;
  //             modphase1.x[] = HUGE;
  //             modphase0.x[] = HUGE;
  //             Point neib = neighborp(1);
  //             if(is_boundary(neighbor(1))){
  //                  foreach_neighbor(1){
  //                     if(point.i==neib.i && point.j==neib.j && point.k==neib.k){
  //                        modphase1.x[] = HUGE;
  //                        modphase0.x[] = HUGE;
  //                        complete1.x[] = 0;
  //                     }
  //                  }
  //             }
  //         }
  //      }
  // }

  foreach_face(){
          complete1.x[] = 0;
          modphase1.x[] = HUGE;
          modphase0.x[] = HUGE;
  }
  // for(int l=0;l<=depth();l++){
  //    foreach_level(l){
  //        foreach_dimension(){
          foreach_face(){
               if((ff4[-1]==0.0 && ff4[]==1.0) || (ff4[-1]==1.0 && ff4[]==0.0)){
                     modphase1.x[]=0.5;
                     modphase0.x[]=0.5;
                     complete1.x[] = 6;
                }else if( (is_phase(ff4[-1]) && !is_phase(ff4[])) || (!is_phase(ff4[-1]) && is_phase(ff4[]))){ // 
                    complete1.x[]=3;
		                 if((is_phase(ff4[-1]) && !is_phase(ff4[])) && orientation(hhh.x[])==0){ // get from height,(1 0)           
                        if(height(hhh.x[])>-1.0 && height(hhh.x[])<0.0){
		                      modphase0.x[]=fabs(height(hhh.x[]));
                          modphase1.x[]=1.0-modphase0.x[];
                          complete1.x[]=2;
                         }
                       
		                 }else if((!is_phase(ff4[-1]) && is_phase(ff4[])) && orientation(hhh.x[])==1){ // 0 1
	                       if(height(hhh.x[])<0.0 && height(hhh.x[]>-1.0)){
                           modphase1.x[]=fabs(height(hhh.x[]));
                           modphase0.x[]=1.0-modphase1.x[];
                           complete1.x[]=2;
                          }
		                 }
                 if(complete1.x[]==3){ // get from middle cell
		                 if((is_phase(ff4[-1]) && !is_phase(ff4[])) && orientation(hhh.x[-1])==0){ // get from height,(1 0)           
                        if(height(hhh.x[-1])>0.0 && height(hhh.x[-1])<1.0){
		                      modphase1.x[]=fabs(height(hhh.x[-1]));
                          modphase0.x[]=1.0-modphase1.x[];
                          complete1.x[]=2;
                         }
                       
		                 }else if((!is_phase(ff4[-1]) && is_phase(ff4[])) && orientation(hhh.x[-1])==1){ // 0 1
	                       if(height(hhh.x[-1])>0.0 && height(hhh.x[-1]<1.0)){
                           modphase0.x[]=fabs(height(hhh.x[-1]));
                           modphase1.x[]=1.0-modphase0.x[];
                           complete1.x[]=2;
                          }
		                 }
                   }
                 if(complete1.x[]==2 && (modphase1.x[]<lim_cut)){
                             modphase1.x[] = lim_cut;
                             modphase0.x[] = 1.0 - modphase1.x[];
                             complete1.x[] = 5;
                    }else if(complete1.x[]==2 && (modphase1.x[]>1.0-lim_cut)){
                             modphase1.x[] = 1.0 - lim_cut;
                             modphase0.x[] = 1.0 - modphase1.x[];
                             complete1.x[] = 5;
                   } 
                  
	           }else{
                  modphase1.x[] = 1.0;
                  modphase0.x[] = 1.0;
                  complete1.x[] = 7;
	          }
        }
    foreach_face (x){
          if(complete1.x[]==3 && modphase1.x[]>HUGE/2.0){
               double c_temp[3][3];
               for(int i0=-1;i0<2;i0++){
                  for(int j0=-1;j0<2;j0++){
                            //coord leftvolume,rightvolume;  //different from linear2.h, rifhtvolume is the right volume of -1; left volume if left volume of 0;
                            double leftvolume,rightvolume;
                            coord a_left,b_left;
                            coord a_right,b_right;
                            coord n_temp;
                            n_temp.x = nn4.x[i0,j0], n_temp.y = nn4.y[i0,j0];

                            coord n_temp0;
                            n_temp0.x = nn4.x[i0-1,j0], n_temp0.y = nn4.y[i0-1,j0];
                            double alpha_temp = alpha4[i0,j0];
                            double alpha_temp0 = alpha4[i0-1,j0];
                            //foreach_dimension(){
                                
                                a_left=(coord){-0.5,-0.5};
                                b_left=(coord){0.5,0.5};
                                
                                    
                                a_right=(coord){-0.5,-0.5};
                                b_right=(coord){0.5,0.5};   
                              
                                b_left.x = 0.0;
                                leftvolume=rectangle_fraction(n_temp,alpha_temp,a_left,b_left)*(b_left.x-a_left.x);        
                                a_right.x = 0.0;
                                rightvolume=rectangle_fraction(n_temp0,alpha_temp0,a_right,b_right)*(b_right.x-a_right.x);  
                                                  
                                c_temp[i0+1][j0+1] = rightvolume + leftvolume;
                          //}// foreach_dimension
                  } //for j0
                } //for i0
                //I am here
                coord n_middle = mycs2(c_temp);//in myc2d.h //mycs3(c_temp);  //get normal from middle volume distribution
                double c_middle = c_temp[0+1][0+1];
                double alpha_middle = plane_alpha(c_middle,n_middle);
                double test=0;
                if(fabs(n_middle.x) > 1e-12){
	           test =  alpha_middle/n_middle.x - (-0.5);
		   // fprintf(fp18,"nnx=%g nny=%g test=%g\n",nn.x,nn.y,test);
                } 
           //     double lim_cut = 1e-13;
                if(test>lim_cut && test<1-lim_cut)
                  {
                        test = test;
                  }else if(test>1.0-lim_cut){
                        test = 1.0 - lim_cut;
                  }else if(test<lim_cut){
                        test = lim_cut;
                  }

                modphase1.x[]=test;
                modphase0.x[]=1.0-modphase1.x[];  
                //complete.x = 4; 
                complete1.x[] = 4;           
   //
          }//complete.x==3
  } //foreach_face
    foreach_face (y){
          if(complete1.y[]==3 && modphase1.y[]>HUGE/2.0){
               double c_temp[3][3];
               for(int i0=-1;i0<2;i0++){
                  for(int j0=-1;j0<2;j0++){
                          //coord leftvolume,rightvolume;  //different from linear2.h, rifhtvolume is the right volume of -1; left volume if left volume of 0;
                          double leftvolume,rightvolume;
                          coord a_left,b_left;
                          coord a_right,b_right;
                              //coord n_temp= interface_normal (point, ff);
                          coord n_temp;
                          n_temp.x = nn4.x[i0,j0], n_temp.y = nn4.y[i0,j0];
                          coord n_temp0;
                          n_temp0.x = nn4.x[i0,j0-1], n_temp0.y = nn4.y[i0,j0-1];
                              //double alpha_temp = plane_alpha (ff[], n_temp);
                          double alpha_temp = alpha4[i0,j0];
                          double alpha_temp0 = alpha4[i0,j0-1];
                              //foreach_dimension(){
                              
                          a_left=(coord){-0.5,-0.5};
                          b_left=(coord){0.5,0.5};
                          
                              
                          a_right=(coord){-0.5,-0.5};
                          b_right=(coord){0.5,0.5};
                        
                          b_left.y = 0.0;
                          leftvolume=rectangle_fraction(n_temp,alpha_temp,a_left,b_left)*(b_left.y-a_left.y);        
                          a_right.y = 0.0;
                          rightvolume=rectangle_fraction(n_temp0,alpha_temp0,a_right,b_right)*(b_right.y-a_right.y);  
                                            
                          c_temp[i0+1][j0+1] = rightvolume + leftvolume;
                            //}// foreach_dimension
                  } //for j0
                } //for i0
                coord n_middle = mycs2(c_temp);//in myc2d.h //mycs3(c_temp);  //get normal from middle volume distribution
                double c_middle = c_temp[0+1][0+1];
                double alpha_middle = plane_alpha(c_middle,n_middle);
                double test=0;;
                if(fabs(n_middle.y) > 1e-12){
	                   test =  alpha_middle/n_middle.y - (-0.5);
		                 // fprintf(fp18,"nnx=%g nny=%g test=%g\n",nn.x,nn.y,test);
                } 
             //   double lim_cut = 1e-13;
                if(test>lim_cut && test<1-lim_cut)
                  {
                        test = test;
                  }else if(test>1.0-lim_cut){
                        test = 1.0 - lim_cut;
                  }else if(test<lim_cut){
                        test = lim_cut;
                  }

                modphase1.y[]=test;
                modphase0.y[]=1.0-modphase1.y[];
               // complete.y = 4; 
                complete1.y[] = 4;               
   //           
          }// complete.y==3
  } //foreach_face (y)       
        //  foreach_dimension(){
        //     if(is_boundary(neighbor(1))){
        //         modphase1.x[]=
        //     }
        //  }


   //set boundary for mod
    foreach_face() {
         
    }
     	  foreach_face(){
                     if(fabs(modphase0.x[])<lim_cut){
                          modphase0.x[]=lim_cut;
                          modphase1.x[]=1.0-modphase0.x[];
                     }else if(fabs(modphase1.x[])<lim_cut){
                          modphase1.x[]=lim_cut;
                          modphase0.x[]=1.0-modphase1.x[];
                     }
                     if(modphase0.x[]>1){
                          modphase0.x[] = 1.0;
                     }
                     if(modphase1.x[]>1){
                          modphase1.x[] = 1.0;
                     }
         }
    foreach_dimension(){
      modphase0.x.restriction = face_rejection2_x ;
      //   modphase0.x.prolongation = mod_prol_x;
         modphase1.x.restriction = face_rejection2_x ;
        // modphase0.x.prolongation = mod_prol_x;
    }
       restriction({modphase0.x, modphase0.y, modphase1.x,modphase1.y});
    
    //  char name101[80];
	  //  sprintf(name101,"mass_record5-pid%d.dat",pid());
	  //   FILE * fp101 = fopen(name101,"w");
    //   foreach_face(){
    //      // if(fabs(modphase1.x[])<=1.0){
    //             fprintf(fp101,"%g %g %g %g %g\n",x,y,z,modphase1.x[],modphase0.x[]);
    //       //}
    //   }
    // fclose(fp101);
    //   MPI_Barrier(MPI_COMM_WORLD);
    //   if(pid()==0){
    //       char command1[150];
    //       sprintf(command1, "LC_ALL=C cat mass_record5-pid*.dat > outfacets/mass_record5-%g",t);
    //       system(command1);

    //       char command7[150];
    //       sprintf(command7, "LC_ALL=C rm -rf mass_record5-pid*.dat");
    //       system(command7);
    //   }
  
  }



///for masstr[] source_pc[] 
//scalar masstr[];
extern scalar masstr;
extern scalar source_pc;
extern scalar source_pc2;
extern scalar vtr;
extern double Tkg,Tkl,Trhog,Trhol,Tcpg,Tcpl,hfg;
extern double source_total,total_area;

extern scalar phase0Tgrad,phase1Tgrad;
//extern scalar source_ff1,source_ff2,ff_old1;
//#if TGRAD_LEON
#if 1
     #include "Tgrad-leon.h"
#endif 


void mass_transfer_rate(){
   //scalar masstr[];
   
   source_total= 0.0;
   total_area=0.0;
  // scalar phase0Tgrad[],phase1Tgrad[];
//#if TGRAD_LEON
#if 1
      Tgrad_leon(ff,phase0Tgrad,phase1Tgrad);
      printf("Trad_leon\n");
// #else
//        T_gradient_dirichlet(phase0Tgrad,phase1Tgrad);
//        printf("Tgrad_classic");
#endif

/*
char name29[80];
sprintf(name29,"gradTwhole%g.dat", t);
FILE * fp29 = fopen(name29,"a");
//fprintf(fp29,"#I J ff Phase0Tgrad Phase1grad temp6\n");
foreach(){
  if(fabs(phase1Tgrad[])>0.0 || fabs(phase0Tgrad[])>0.0){
     double temp6;
     temp6 = phase0Tgrad[] - phase1Tgrad[]; 
     fprintf(fp29,"%d %d %g %g %g %g\n",point.i,point.j,ff[],phase0Tgrad[],phase1Tgrad[],temp6);
  }
}
fclose(fp29);
*/
//   if((globali%outstep)==0){
        // char name33[80];
        // sprintf(name33,"mass_record3-pid%d.dat",pid());
        // FILE * fp33 = fopen(name33,"w");
        // foreach(){
        //     fprintf(fp33,"%g %g %g %g %g\n",x,y,T[],phase0Tgrad[],phase1Tgrad[]);
        // }
        // fclose(fp33);
        // MPI_Barrier(MPI_COMM_WORLD);
        // if(pid()==0){
        //                 char command1[150];
        //                 sprintf(command1, "LC_ALL=C cat mass_record3-pid*.dat > outfacets/mass_record3-%g",t);
        //                 system(command1);

        //                 char command7[150];
        //                 sprintf(command7, "LC_ALL=C rm -rf mass_record3-pid*.dat");
        //                 system(command7);
        //     }
//   }


   foreach(reduction(+:source_total) reduction(+:total_area)){
         masstr[]=0.0;
         source_pc[] = 0.0;
        // source_ff1[] = 0.0;
     // fixed me; when ff[]=1.0 and ff[neighbour]=0.0, phase0Tgrad is failed
     // fixed me: boundary for source_pc and masstr 
        if(ff[]<1.0-EPS && ff[]>EPS){
            masstr[] = Tkg/hfg*phase0Tgrad[] - Tkl/hfg*phase1Tgrad[];
		 
	 }else{
	      masstr[]=0.0;
	 }
     // caculate surface area (non - dimension!!!!!!!!!!!!!!)
       double area = 0.;
        coord pp = {0.,0.};
	//  foreach (reduction(+:area))
	    if (ff[] > EPS && ff[] < 1. - EPS) {
	      coord n = interface_normal (point, ff);
          
	      double alpha = plane_alpha (ff[], n);
               area = line_length_center(n, alpha, &pp);//plane_area_center (n, alpha, &p);
	      //area += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
	    }

        //s_v(i,j,k) = mdot(i,j,k)*inv_drho*cell_area(i,j,k)/dh //area is a ratio, not a dimension number
                 source_pc[] = masstr[]*(1.0/Trhog - 1.0/Trhol)*area/Delta;
        #if AXI
                  source_pc2[] = source_pc[]*max(1e-20,(y+pp.y*Delta));
                //  source_pc2[] = 0.0; //source_pc[]*(y+pp.y*Delta); //20230603
                // printf("source_pc9999\n");
        #endif
        source_total = source_total + source_pc[];
        total_area = total_area + area;

       // source_ff1[] = 1; //cells contain sources;
   }

 //if((globali%outstep)==0){
  /*
if(out_flag){
  get_T_from_Tlff_Tgff3(1,ff);
	char name59[80];
	sprintf(name59,"outfacets/mass_record1-pid%d.dat", pid());
	FILE * fp59 = fopen(name59,"w");
       // fprintf(fp59,"#i j masstr source_pc\n");
        foreach(){
          
            // if(ff[]<1.0-EPS && ff[]>EPS){
         //   if(ff[]<1.0-CCCCMIN && ff[]>CCCCMIN){
                coord n = interface_normal (point, ff), p;
	        double alpha = plane_alpha (ff[], n);
                double area = plane_area_center (n, alpha, &p);
                double cha =phase0Tgrad[] -phase1Tgrad[];
                if(level == level_interface && (fabs(z-L0/2.0)<=Delta/2.0)){
                  	 //fprintf(fp59,"%g %g %g %g %g %g %g %g %g %g %d %g\n",x,y,ff[],masstr[],source_pc[],T[],area,phase0Tgrad[],phase1Tgrad[],cha,level,Delta);
                     fprintf(fp59,"%g %g\n",x,y);
                 }
      //       }
             
             //if(fabs(source_pc[])>EPS)
             // fprintf(fp59,"%g %g %g %g %g %d %g\n",x,y,ff[],source_pc[],T[],level,Delta);
        }
        fclose(fp59);
  //       char name60[80];
	// sprintf(name60,"outfacets/source_area_record.dat");
	// FILE * fp60 = fopen(name60,"a");
  //       fprintf(fp60,"%g %g %g\n",t,source_total,total_area);
  //       fclose(fp60);


            char name_temp[80];
    sprintf(name_temp,"outfacets/mass_record1-pid%d.dat", pid());
    mass_record_out(name_temp);

            char command1[150];
          sprintf(command1, "LC_ALL=C cat outfacets/mass_record1-pid*.dat > outfacets/mass_record1-%g",t);
     //     system(command1);

          char command7[150];
          sprintf(command7, "LC_ALL=C rm -rf  outfacets/mass_record1-pid*.dat");
     //     system(command7);
}
    
*/
     foreach(){
         vtr[]=0.0;
         if(fabs(masstr[])>0.0){
            vtr[]=masstr[]/Trhol;
         }
     }


 
}



void LevelSetShift2VOFChange(double dt){
    
     scalar overshoot[];
     foreach(){
          deltac[] = 0.0;
     }
int flag_overshoot=0;
     foreach(reduction(max:flag_overshoot)){
              overshoot[]=0.0;
          // if(ff[]>EPS && ff[]<1.0-EPS){
	     if(ff[]>1e-6 && (ff[]<1-1e-6) && (css_test3_n[]>0.0)){
              coord n= interface_normal (point, ff);
	            double alpha_old = plane_alpha (ff[], n);
              double magn=sqrt(n.x*n.x + n.y*n.y);// + n.z*n.z); //2D-3D
              double magn1=fabs(n.x)+fabs(n.y);//+fabs(n.z);
              double nshift = masstr[]*dt/Trhol/Delta; 


	      double dalpha=-nshift*magn; //
              double alpha_new = alpha_old;
              //if(alpha_old + dalpha > 1.0){  //this is wrong, because -0.5<alpha<0.5
              if(alpha_old + dalpha > 0.5){
                 alpha_new = 0.5;
                 overshoot[] = alpha_old + dalpha - 0.5;
                 deltac[] = deltac[] + (1.0 - ff[]);
              }else if(alpha_old + dalpha < -0.5){
	               alpha_new = - 0.5;
                 overshoot[] = alpha_old + dalpha - (-0.5);
                 deltac[] = deltac[] - ff[];
              }else{
		             alpha_new = alpha_old + dalpha;
                 double ffnew = plane_volume (n,alpha_new);
                 deltac[] = deltac[] + (ffnew - ff[]);
              }
              
  //            fprintf(fp37,"i=%d j=%d f=%g deltac=%g dalpha=%g alpha_old=%g  n.x=%g n.y=%g abs2(n)=%g abs1(n)=%g overshoot=%g \n",point.i,point.j,ff[],deltac[],dalpha,alpha_old,n.x,n.y,magn,magn1,overshoot[]);
            //  if(fabs(overshoot[]*magn)>0.3){ //overshoot*magn should be smaller than 0.1
	      //  if(fabs(overshoot[])>0.1){
             //        flag_overshoot=1;
             //   }
           }
     }


 //    fclose(fp37);
   //  fclose(fp40)
     if(flag_overshoot == 1){
          fprintf(stderr,"overshoot*magn >= 0.1 \n");
          exit(1);
      }


     scalar intmask[];
     foreach(){
	      intmask[] = 3;
        if(ff[]<1.0-EPS && ff[]>EPS){
 		       intmask[] = 0;
        }
     }
     foreach(){
	      bool is1= false;
        foreach_dimension(){
	      if(((intmask[1]==0) || (intmask[-1]==0)) && level==level_interface){
	               is1 = true;
           }
        }
        if (is1 && intmask[]!=0){
           intmask[] = ff[]>0.5-EPS3 ? 1 : -1 ;
        }
        
     }

/*
char name62[80];
sprintf(name62,"intmask%g.dat", t);
FILE * fp62 = fopen(name62,"a");
fprintf(fp62,"#i j intmask\n");
foreach(){
 double rr=sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5));
 if(rr<0.01+0.2 && rr>0.2-0.01 )
    fprintf(fp62,"%d %d %g\n",point.i,point.j,intmask[]);
  
}
fclose(fp62);     
*/
     //handle overshoot; because it's difficut to calculate neighbour's normal; so this version is opposite to paris that we first search overshoot, then judge which cell accept this overshoot
     double normal_limit; // this should be think more

   normal_limit = 1e-3;
   /*
     if(stefan_flat!=1){
          normal_limit = 1e-3;
     }else{
          normal_limit = 1e-1;
     }
   */

     vector position2[];
     scalar deltac_temp[];
//char name42[80];
//sprintf(name42,"drainwarning-%g.dat", t);
//FILE * fp42 = fopen(name42,"a");

//char name46[80];
//sprintf(name46,"overshootbig-%g.dat", t);
//FILE * fp46 = fopen(name46,"a");

     foreach(){
       
        deltac_temp[] = 0.0;
        if(level==level_interface){
        position2.x[]=0;
        position2.y[]=0;
       // position2.z[]=0;
        
	//if(fabs(overshoot[])>0.0){
  	if(fabs(overshoot[])>0.0 && (css_test3_n[]>0.0)){
  //            fprintf(fp46,"i=%d j=%d overshoot=%g ",point.i,point.j,overshoot[]);
              coord accept_possition;
//              accept_possition.x = 0;
//              accept_possition.y = 0;
              coord n = interface_normal (point, ff);
/*
	      if(fabs(n.x)>normal_limit){
		  accept_possition.x = sign(n.x);
              }
              if(fabs(n.y)>normal_limit){
		  accept_possition.y = sign(n.x);
              }
*/           
/* 
              foreach_dimension(){
		if(fabs(n.x)>normal_limit){
		    accept_possition.x = sign(n.x);
                 }
              }
*/
//              int xx = accept_possition.x;
//              int yy = round(accept_possition.y);
              int xx=0,yy=0;//,zz=0;
              if(fabs(n.x)>normal_limit){
		              if(overshoot[]>0.0){
                    xx = sign(n.x); 
                  }else{
                    xx = -sign(n.x);
                  }
               }
              if(fabs(n.y)>normal_limit){
                  if(overshoot[]>0.0){  
                        yy = sign(n.y);        
                  }else{
                        yy = -sign(n.y);
                  }
               }
            //    if(fabs(n.z)>normal_limit){
            //       if(overshoot[]>0.0){  
            //             zz = sign(n.z);        
            //       }else{
            //             zz = -sign(n.z);
            //       }
            //    }
              if(xx==0 && yy==0){// && zz==0){
                   fprintf (stderr, "WARNING: xx=0 && yy=0 && zz==0\n");
               }
 // for boiling case below should be take care!!!!!!!!!!!!!!!!!!!!!!!    
 //mpi point.i take care
      // if(point.i+xx>BGHOSTS2 && point.j+yy>BGHOSTS2 && point.i+xx<N+BGHOSTS2+1 && point.j+yy<N+BGHOSTS2+1){
      //   if(!is_boundary_box2(neighbor(xx,yy,zz))){
     // if(!is_boundary_box2(neighbor(xx,yy,zz)) && (css_test3_n[])>0.0){
     if(!is_boundary_box2(neighbor(xx,yy)) && (css_test3_n[])>0.0){
            //  if(fabs(intmask[xx,yy,zz]-3.0)>0.5){
                 if(fabs(intmask[xx,yy]-3.0)>0.5){
		  // double magn=sqrt(n.x*n.x+n.y*n.y);
                   double alpha_old = plane_alpha (ff[xx,yy], n); //using normal of interface and volume of accept cell
                   //double alpha_new = alpha_old + overshoot[];
                   double alpha_new = (alpha_old + overshoot[]);
                   //if(alpha_new>1.0){ //this is wrong, because -0.5<alpha<0.5
                   if(alpha_new>0.5){
                      //  double temp1 = alpha_new - alpha_old;
                        fprintf (stderr, "WARNING: Overfilling a cell already accommodating overfilling t=%g alp_new=%g ff=%g \n",t,alpha_new,ff[]);
                     //   fprintf (fp42, "WARNING: Overfilling a cell already accommodating overfilling t=%g alp_new=%g ff=%g over=%g\n",t,alpha_new,ff[],temp1);
                        
                        alpha_new = 0.5;
                        int flag3=0;
                        double shoot_cut=1e-12;
                       // for 3d - I don't have the following adaptation; just comment
                      //   if(fabs(alpha_old + overshoot[]-0.5)>shoot_cut){
                      //      int xx2=xx;
                      //      int yy2=yy;
                      //      int zz2=zz; //for 3D 
                      //      //if(fabs(n.x)>=fabs(n.y) && yy!=0 && xx!=0 && ff[xx,0]<=0.0){
                      //      if(fabs(n.x)>=fabs(n.y) && yy!=0 && xx!=0 && ff[xx,0]<ff[xx,yy]){
			                //         	yy2 = 0;
                      //           flag3=1;
                      //      //}else if(fabs(n.x)<fabs(n.y) && xx!=0 && yy!=0 && ff[0,yy]<=0.0){
                      //      }else if(fabs(n.x)<fabs(n.y) && xx!=0 && yy!=0 && ff[0,yy]<ff[xx,yy]){
                      //           xx2 = 0;
                      //           flag3 = 1;
                      //      }

                      //      if(flag3 == 1){
                      //        double alpha_old2 = plane_alpha (ff[xx2,yy2], n);
                      //        double alpha_new2 = alpha_old2 + overshoot[];
                      //        if(alpha_new2<0.5 && alpha_new>-0.5){
                      //           xx = xx2;
                      //           yy = yy2;
                      //           alpha_old = alpha_old2;
                      //           alpha_new = alpha_new2;
                      //        }
                      //      }else{
                      //          fprintf(stderr,"could not do someting\n");
                      //  //        fprintf (fp42,"could not do someting\n");
                      //      }
                             
                      //    }
                        
                       
                   }else if(alpha_new < -0.5){
                        //double temp1 = alpha_new - alpha_old;
                        fprintf (stderr, "WARNING: Draining a cell already accommodating overdrainage! t=%g alp_new=%g ff=%g \n", t,alpha_new,ff[]);
                         
                    //    fprintf (fp42, "WARNING: Draining a cell already accommodating overdrainage! t=%g alp_new=%g ff=%g over=%g\n", t,alpha_new,ff[],temp1);
                        alpha_new = -0.5;
                    //add by xiangbin 2022-0801
                        int flag3=0;
                        double shoot_cut=1e-12;
                         // for 3d - I don't have the following adaptation; just comment
                    //     if(fabs(alpha_old + overshoot[]-(-0.5))>shoot_cut){
                    //        int xx2=xx;
                    //        int yy2=yy; 
                    //        //if(fabs(n.x)>=fabs(n.y) && yy!=0 && xx!=0 && ff[xx,0]>=1.0){
                    //        if(fabs(n.x)>=fabs(n.y) && yy!=0 && xx!=0 && ff[xx,0]>ff[xx,yy]){
				            //             yy2 = 0;
                    //             flag3=1;
                    //        //}else if(fabs(n.x)<fabs(n.y) && xx!=0 && yy!=0 && ff[0,yy]>=1.0){
                    //        }else if(fabs(n.x)<fabs(n.y) && xx!=0 && yy!=0 && ff[0,yy]>ff[xx,yy]){
                    //             xx2 = 0;
                    //             flag3 = 1;
                    //        }

                    //        if(flag3 == 1){
                    //          double alpha_old2 = plane_alpha (ff[xx2,yy2], n);
                    //          double alpha_new2 = alpha_old2 + overshoot[];
                    //          if(alpha_new2<0.5 && alpha_new>-0.5){
                    //             xx = xx2;
                    //             yy = yy2;
                    //             alpha_old = alpha_old2;
                    //             alpha_new = alpha_new2;
                    //          }
                    //        }else{
                    //            fprintf(stderr,"could not do someting\n");
                    //  //          fprintf (fp42,"could not do someting\n");
                    //        }
                             
                    //      }
                   }
              
           // for boiling case above should be take care!!!!!!!!!!!!!!!!!!!!!!!
             double ffnew = plane_volume (n,alpha_new);
              //deltac[xx,yy] = deltac[xx,yy] + cnew - ff[xx,yy];
              deltac_temp[] =  ffnew - ff[xx,yy];
              position2.x[]=xx;
              position2.y[]=yy;
             // position2.z[]=zz;
            }

           }
        }
       }//level == level_interface
     }
//fclose(fp42);
//fclose(fp46);

//char name38[80];
//sprintf(name38,"i0j0i1j1-%g.dat", t);
//FILE * fp38 = fopen(name38,"a");
scalar flagi[];
foreach(){
   flagi[]=0;
}
foreach(){
  //if(fabs(intmask[]-3.0)>0.5){
   if(fabs(intmask[]-3.0)>0.5 && css_test3_n[]>0.0){
    for(int i0=-1; i0<=1; i0=i0+1){
       for(int j0=-1; j0<=1; j0=j0+1){            
        //  for(int k0=-1; k0<=1; k0=k0+1){
        //if neighbor has overshoot, judge if that overshoot belongs toyou
                  if(!(i0==0 && j0==0) ){
                      if(fabs(deltac_temp[i0,j0])>0.0){
                                int i1=-round(position2.x[i0,j0]);
                                int j1=-round(position2.y[i0,j0]);   
                                if(i0==i1  && j0==j1 && (css_test3_n[i0,j0]>0.0)) {	
                                      deltac[] = deltac[] + deltac_temp[i0,j0];
                                              // flagi[i0,j0]=flagi[i0,j0]+1;
                                              flagi[] = flagi[] + 1;
                            //                  fprintf(fp38,"i=%d j=%d ff=%g i0=%d j0=%d i1=%d j1=%d pos2(i0,j0).x=%g pos2(i0,j0).y=%g deltac=%g deltac_temp(i0,j0)=%g \n",point.i,point.j,ff[],i0,j0,i1,j1,position2.x[i0,j0],position2.y[i0,j0],deltac[],deltac_temp[i0,j0]);
                                }
                      }
                 }
         // } //k0
       } //j0
    } //i0
  }
}


}


void  mov_interface_dc(){

double small_fraction_limit_T = 1e-12; //EPS;
    // foreach(){
    //     if(ff[]<small_fraction_limit_T && ff[]>0.0){ //this is just for evaporate, if cooling add >1e-6 situation
    //        Tlff[]=0.0;
    //        Tgff[]=Tgff[]/(1.0-ff[])*(1.0-0.0);
    //        ff[]=0.0;
    //     }else if(ff[]>1-small_fraction_limit_T && ff[]<1.0){
    //        Tlff[]=Tlff[]/ff[]*1.0;
    //        Tgff[]=0.0;
    //        ff[]=1.0;
    //     }
    // }

scalar f_old[];
foreach(){
      if(css_test3_n[]<1.0){ //css_test3[]>0.0
          T_solid[]=T[];
      }
   }



foreach(){
//move dc + clamp(f[],0.0,1,0)
        f_old[] = ff[];
        ff[] = ff[] + deltac[];
	    
   }
for(scalar s in {T}){
      foreach(){
           double val_tot=0.0;
           double wei=0.0;
           double val1=0.0;
           if(css_test3_n[]>0.0){
                    // if(ff[] > 0.0 && f_old[] <= 0.0 ){ //
                    //     val_tot = val_tot + Tsat00*ff[];
                    //     wei = wei + ff[];
                    //     val_tot = val_tot + T[]*(1.0-ff[]);
                    //     wei = wei + (1.0-ff[]);
                    // }else if(1.0-ff[] > 0.0 && 1.0-f_old[] <= 0.0){
                    //     val_tot = val_tot + Tsat00*(1.0-ff[]);
                    //     wei = wei + (1.0-ff[]);
                    //     val_tot = val_tot + T[]*ff[];
                    //     wei = wei + ff[];
                    // }else{
                    //     val_tot = s[];
                    //     wei = 1.0;
                    // }
                    if(ff[] > 0.5 && f_old[] <= 0.5 ){ // 20230605
                        val_tot = val_tot + Tsat00*ff[];
                        wei = wei + ff[];
                        val_tot = val_tot + T[]*(1.0-ff[]);
                        wei = wei + (1.0-ff[]);
                    }else if(1.0-ff[] > 0.5 && 1.0-f_old[] <= 0.5){
                        val_tot = val_tot + Tsat00*(1.0-ff[]);
                        wei = wei + (1.0-ff[]);
                        val_tot = val_tot + T[]*ff[];
                        wei = wei + ff[];
                    }else{
                        val_tot = s[];
                        wei = 1.0;
                    }
                    val1 = val_tot/wei;
           }
           if(css_test3_n[]<1.0){ //css_test3[]>0.0
                s[] = val1*css_test3_n[] + T_solid[]*(1.0-css_test3_n[]);
            }else{
                s[] = val1;
            }
	    }

}
boundary({T}); //fiexed me: need to consider boundary condition for T embed;

// foreach(){
// 	    for(scalar s in {Tlff,Tgff}){
// 		if(s.inverse){
// 		     s[]=s[]*(1.0-ff[]);
// 		 }else{
// 	             s[]=s[]*ff[];
// 		 }
// 	    }
//    }

}


void ulf_ugf_function(scalar topo_mask, scalar topo_mask_g){
    face vector facemask[];
    foreach_face(){
        facemask.x[]=1;
        //if(is_boundary_box2(point) || is_boundary_box2(neighbor(-1))){   
        //     facemask.x[]=0;
        //}
        if(is_boundary_box2(cell) || is_boundary_box2(neighbor(-1))){   
             //facemask.x[]=0;
             facemask.x[]=1;
        }
     }

     foreach_face(){
	       ulf.x[] = 0.0;
         ugf.x[] = 0.0;
     }
         
     
     //////////////////////////////////////////////////////////
     // ???assume that interface cell will not be coarsen or refine
     //////////////////////////////////////////////////////////
     foreach_face(){
         //vofnonmodule.f90 165-193
       if(facemask.x[]==1){
          if (level==level_interface){
              int correct1 =0;
              if((topo_mask[-1]<1 && topo_mask[]<1) && (topo_mask[]==0 || topo_mask[-1]==0)  ){
                      ulf.x[] = uf.x[] + usf.x[];
                }else{
                      ulf.x[] = uf.x[];
                }
          }else{
              ulf.x[] = uf.x[];
          }

        }
        #if EMBED
         if(fabs(fs.x[])<=0){
             ulf.x[] = 0.0; 
             usf.x[] = 0.0;
             uf.x[] = 0.;
         }
        #endif
     }

    

    foreach_face(){
         //vofnonmodule.f90 165-193
       if(facemask.x[]==1){
          if (level==level_interface){
              int correct1 =0;
              if((topo_mask_g[-1]<1 && topo_mask_g[]<1) && (topo_mask_g[]==0 || topo_mask_g[-1]==0)  ){
              // if(fabs(topo_mask[])<=2){ 
                      ugf.x[] = uf.x[] + usfg.x[];
                }else{
                      ugf.x[] = uf.x[];
                }
          }else{
              ugf.x[] = uf.x[];
          }

        }
        #if EMBED
        if(fabs(fs.x[])<=0){
            ugf.x[] = 0.0; 
            usfg.x[] = 0.0;
            uf.x[] = 0.;
        }
        #endif
         //if(fabs(topo_mask[])<3)
 		          //fprintf(fp80,"%g %g %g %g %g\n",x,y,ulf.x[],uf.x[],usf.x[]);
     }
    // // // if((globali%outstep)==0){
    // // //       char name80[80];
    // // //       sprintf(name80,"outfacets/ulf-check-%g.dat",t);
    // // //       FILE * fp80 = fopen(name80,"a");
    // // //     foreach_face(x){
    // // //         if(fabs(usf.x[])>0.0)
    // // //               fprintf(fp80,"%g %g %g %g %g\n",x,y,ulf.x[],uf.x[],usf.x[]);      
    // // //     }
    // // //     fclose(fp80);
    // // // }
}

static void embed_fraction_refine2 (Point point, scalar css_test)
//void embed_fraction_refine2 (Point point, scalar css_test)
{
  double cc = css_test[];

  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (cc <= 0. || cc >= 1.) {
    foreach_child()
      css_test[] = cc;
  }
  else {

    /**
    If the cell contains the embedded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, css_test, fss_test);
    double alpha = plane_alpha (cc, n);
      
    foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      foreach_dimension()
	nc.x = child.x*n.x;
      css_test[] = rectangle_fraction (nc, alpha, a, b);
    }
  }
}

foreach_dimension()
static void embed_face_fraction_refine2_x (Point point, scalar s)
{
 // printf("YOUYONG\n");
  vector fss_test = s.v;

  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (css_test[] <= 0. || css_test[] >= 1.) {

    /**
    We need to make sure that the fine cells face fractions match
    those of their neighbours. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	fine(fss_test.x,1,j,k) = css_test[];
    for (int i = 0; i <= 1; i++)
      if (!is_refined(neighbor(2*i-1)) && neighbor(2*i-1).neighbors &&
	  (is_local(cell) || is_local(neighbor(2*i-1))))
	for (int j = 0; j <= 1; j++)
	  for (int k = 0; k <= 1; k++)
	    fine(fss_test.x,2*i,j,k) = fss_test.x[i];
  }
  else {

    /**
    If the cell contains the embedded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, css_test, fss_test);
    double alpha = plane_alpha (css_test[], n);
      
    /**
    We need to reconstruct the face fractions *fss_test* for the fine cells.
    
    For the fine face fractions contained within the coarse cell,
    we compute the intersections directly using the VOF
    reconstruction. */

#if dimension == 2

    /**
    In 2D, we obtain the face fractions by taking into
    account the orientation of the normal. */

    if (2.*fabs(alpha) < fabs(n.y)) {
      double yc = alpha/n.y;
      int i = yc > 0.;
      fine(fss_test.x,1,1 - i) = n.y < 0. ? 1. - i : i;
      fine(fss_test.x,1,i) = n.y < 0. ? i - 2.*yc : 1. - i + 2.*yc;
    }
    else
      fine(fss_test.x,1,0) = fine(fss_test.x,1,1) = alpha > 0.;

#else // dimension == 3

    /**
    in 3D, we use the 2D projection of the reconstruction. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++) {
	if (!fine(css_test,0,j,k) || !fine(css_test,1,j,k))
	  fine(fss_test.x,1,j,k) = 0.;
	else {
	  static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
	  coord nc;
	  nc.x = 0., nc.y = (2.*j - 1.)*n.y, nc.z = (2.*k - 1.)*n.z;
	  fine(fss_test.x,1,j,k) = rectangle_fraction (nc, alpha, a, b);
	}
      }

#endif // dimension == 3
    
    /**
    For the fine face fractions coincident with the faces of the
    coarse cell, we obtain the intersection position from the
    coarse cell face fraction. */

    for (int i = 0; i <= 1; i++)
      if (neighbor(2*i-1).neighbors &&
	  (is_local(cell) || is_local(neighbor(2*i-1)))) {
	if (!is_refined(neighbor(2*i-1))) {
	  if (fss_test.x[i] <= 0. || fss_test.x[i] >= 1.)
	    for (int j = 0; j <= 1; j++)
	      for (int k = 0; k <= 1; k++)
		fine(fss_test.x,2*i,j,k) = fss_test.x[i];
	  else {
#if dimension == 2
	  
	    /**
	    In 2D the orientation is obtained by looking at the values
	    of face fractions in the transverse direction. */
	  
	    double a = fss_test.y[0,1] <= 0. || fss_test.y[2*i-1,1] <= 0. ||
	      fss_test.y[] >= 1. || fss_test.y[2*i-1] >= 1.;
	    if ((2.*a - 1)*(fss_test.x[i] - 0.5) > 0.) {
	      fine(fss_test.x,2*i,0) = a;
	      fine(fss_test.x,2*i,1) = 2.*fss_test.x[i] - a;
	    }
	    else {
	      fine(fss_test.x,2*i,0) = 2.*fss_test.x[i] + a - 1.;
	      fine(fss_test.x,2*i,1) = 1. - a;
	    }

#else  // dimension == 3

	    /**
	    In 3D we reconstruct the face fraction from the projection
	    of the cell interface reconstruction, as above. */
	  
	    for (int j = 0; j <= 1; j++)
	      for (int k = 0; k <= 1; k++) {
		static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
		coord nc;
		nc.x = 0., nc.y = (2.*j - 1.)*n.y, nc.z = (2.*k - 1.)*n.z;
		fine(fss_test.x,2*i,j,k) =
		  rectangle_fraction (nc, alpha - n.x*(2.*i - 1.)/2., a, b);
	      }

#endif // dimension == 3
	  }
	}

	/**
	The face fractions of empty children cells must be zero. */
	
	for (int j = 0; j <= 1; j++)
	#if dimension > 2
	  for (int k = 0; k <= 1; k++)
	#endif
	    if (fine(fss_test.x,2*i,j,k) && !fine(css_test,i,j,k))
	      fine(fss_test.x,2*i,j,k) = 0.;
      }
  }
}

//refine for css_test2
static void embed_fraction_refine3 (Point point, scalar css_test2)
//void embed_fraction_refine3 (Point point, scalar css_test2) //2
{
  double cc = css_test2[];

  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (cc <= 0. || cc >= 1.) {
    foreach_child()
      css_test2[] = cc;
  }
  else {

    /**
    If the cell contains the embedded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, css_test2, fss_test2);
    double alpha = plane_alpha (cc, n);
      
    foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      foreach_dimension()
	nc.x = child.x*n.x;
      css_test2[] = rectangle_fraction (nc, alpha, a, b);
    }
  }
}

foreach_dimension()
static void embed_face_fraction_refine3_x (Point point, scalar s)
{
  vector fss_test2 = s.v;

  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (css_test2[] <= 0. || css_test2[] >= 1.) {

    /**
    We need to make sure that the fine cells face fractions match
    those of their neighbours. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	fine(fss_test2.x,1,j,k) = css_test2[];
    for (int i = 0; i <= 1; i++)
      if (!is_refined(neighbor(2*i-1)) && neighbor(2*i-1).neighbors &&
	  (is_local(cell) || is_local(neighbor(2*i-1))))
	for (int j = 0; j <= 1; j++)
	  for (int k = 0; k <= 1; k++)
	    fine(fss_test2.x,2*i,j,k) = fss_test2.x[i];
  }
  else {

    /**
    If the cell contains the embedded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, css_test2, fss_test2);
    double alpha = plane_alpha (css_test2[], n);
      
    /**
    We need to reconstruct the face fractions *fss_test2* for the fine cells.
    
    For the fine face fractions contained within the coarse cell,
    we compute the intersections directly using the VOF
    reconstruction. */

#if dimension == 2

    /**
    In 2D, we obtain the face fractions by taking into
    account the orientation of the normal. */

    if (2.*fabs(alpha) < fabs(n.y)) {
      double yc = alpha/n.y;
      int i = yc > 0.;
      fine(fss_test2.x,1,1 - i) = n.y < 0. ? 1. - i : i;
      fine(fss_test2.x,1,i) = n.y < 0. ? i - 2.*yc : 1. - i + 2.*yc;
    }
    else
      fine(fss_test2.x,1,0) = fine(fss_test2.x,1,1) = alpha > 0.;

#else // dimension == 3

    /**
    in 3D, we use the 2D projection of the reconstruction. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	if (!fine(css_test2,0,j,k) || !fine(css_test2,1,j,k))
	  fine(fss_test2.x,1,j,k) = 0.;
	else {
	  static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
	  coord nc;
	  nc.x = 0., nc.y = (2.*j - 1.)*n.y, nc.z = (2.*k - 1.)*n.z;
	  fine(fss_test2.x,1,j,k) = rectangle_fraction (nc, alpha, a, b);
	}

#endif // dimension == 3
    
    /**
    For the fine face fractions coincident with the faces of the
    coarse cell, we obtain the intersection position from the
    coarse cell face fraction. */

    for (int i = 0; i <= 1; i++)
      if (neighbor(2*i-1).neighbors &&
	  (is_local(cell) || is_local(neighbor(2*i-1)))) {
	if (!is_refined(neighbor(2*i-1))) {
	  if (fss_test2.x[i] <= 0. || fss_test2.x[i] >= 1.)
	    for (int j = 0; j <= 1; j++)
	      for (int k = 0; k <= 1; k++)
		fine(fss_test2.x,2*i,j,k) = fss_test2.x[i];
	  else {
#if dimension == 2
	  
	    /**
	    In 2D the orientation is obtained by looking at the values
	    of face fractions in the transverse direction. */
	  
	    double a = fss_test2.y[0,1] <= 0. || fss_test2.y[2*i-1,1] <= 0. ||
	      fss_test2.y[] >= 1. || fss_test2.y[2*i-1] >= 1.;
	    if ((2.*a - 1)*(fss_test2.x[i] - 0.5) > 0.) {
	      fine(fss_test2.x,2*i,0) = a;
	      fine(fss_test2.x,2*i,1) = 2.*fss_test2.x[i] - a;
	    }
	    else {
	      fine(fss_test2.x,2*i,0) = 2.*fss_test2.x[i] + a - 1.;
	      fine(fss_test2.x,2*i,1) = 1. - a;
	    }

#else  // dimension == 3

	    /**
	    In 3D we reconstruct the face fraction from the projection
	    of the cell interface reconstruction, as above. */
	  
	    for (int j = 0; j <= 1; j++)
	      for (int k = 0; k <= 1; k++) {
		static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
		coord nc;
		nc.x = 0., nc.y = (2.*j - 1.)*n.y, nc.z = (2.*k - 1.)*n.z;
		fine(fss_test2.x,2*i,j,k) =
		  rectangle_fraction (nc, alpha - n.x*(2.*i - 1.)/2., a, b);
	      }

#endif // dimension == 3
	  }
	}

	/**
	The face fractions of empty children cells must be zero. */
	
	for (int j = 0; j <= 1; j++)
	#if dimension > 2
	  for (int k = 0; k <= 1; k++)
	#endif
	    if (fine(fss_test2.x,2*i,j,k) && !fine(css_test2,i,j,k))
	      fine(fss_test2.x,2*i,j,k) = 0.;
      }
  }
}





// void restriction_Tlff (Point point, scalar s)
// {
//   double sum = 0.;
//   double weight= 0.;
//   foreach_child(){
//     if(ff[]>0.5){
//         sum += s[];
//         weight += ff[];
//     }    
//   }
//   // s[] = sum/(1 << dimension)/(cm[] + 1e-30);
//   if(weight>1e-12){
//     s[] = sum/weight;
//   }
// }

// void restriction_Tgff (Point point, scalar s)
// {
//   double sum = 0.;
//   double weight= 0.;
//   foreach_child(){
//     if(1.0-ff[]>0.5){
//         sum += s[];
//         weight += 1.0-ff[];
//     }    
//   }
//   // s[] = sum/(1 << dimension)/(cm[] + 1e-30);
//   if(weight>1e-12){
//     s[] = sum/weight;
//   }
// }