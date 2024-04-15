#include "heights.h"
#include "fractions.h"

#include "./linear2-tree-2.h"
#include "./vof2front-advanced-copy.h"

// extern face vector modphase0,modphase1;
// extern vector hhh;
extern int globali,outstep,level_interface;
extern scalar css_test3_n;
extern scalar deltac;
extern scalar T_solid;

extern scalar css_test2;
extern face vector fss_test2;

extern scalar css_test;
extern face vector fss_test;
extern bool energy_advecting_flag;
extern bool flag_topos_advect_uf;
extern scalar intersect_true;
extern bool flag_average_source;

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
  // vector hhh;
  face vector modphase1;
  face vector modphase0;
  // scalar Tl;
};

void get_modphase01_2(struct Threephases q){


//check mod//initial
// foreach(){
//     foreach_dimension(){
//       smallmodl.x[]=-10.0;
//       bigmodl.x[]=-10.0;
//       smallmodg.x[]=-10.0;
//       bigmodg.x[]=-10.0;

//     }
// }
///////////////

    
    scalar ff4=q.ff4;
    face vector modphase1=q.modphase1;
    face vector modphase0=q.modphase0;
    // scalar Tl=q.Tl;
    vector hhh[];
    heights(ff4,hhh);
    vector nn4[];
    scalar alpha4[];
    // vector hhh=q.hhh;
    // foreach(){
    //    ff4[]=ff[];
    // }   
    // ff4.restriction = restriction_volume_average;
    // ff4.restriction = restriction_volume_average;
    ff4.refine = ff4.prolongation = fraction_refine;
    
    boundary({hhh});
    
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
          if((complete1.x[]==3 && modphase1.x[]>HUGE/2.0) || complete1.x[] == 5){
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
                              double n_temp0_abs=0;
                            n_temp0_abs = sqrt(n_temp0.x*n_temp0.x+n_temp0.y*n_temp0.y);
                            if(fabs(n_temp0_abs)>0){
                                b_left.x = 0.0;
                                leftvolume=rectangle_fraction(n_temp,alpha_temp,a_left,b_left)*(b_left.x-a_left.x);        
                            }else{
                                if(ff4[i0-1,j0]>=1){
                                    leftvolume=0.5;
                                }else if(ff4[i0-1,j0]<=0){
                                    leftvolume=0;
                                }else{
                                    leftvolume=ff4[i0-1,j0]/2.0;
                                }
                            }
                            double n_temp_abs=0;
                            n_temp_abs = sqrt(n_temp.x*n_temp.x+n_temp.y*n_temp.y);
                            if(fabs(n_temp_abs)>0){
                                a_right.x = 0.0;
                                rightvolume=rectangle_fraction(n_temp0,alpha_temp0,a_right,b_right)*(b_right.x-a_right.x);  
                            }else{
                               if(ff4[i0,j0]>=1){
                                    leftvolume=0.5;
                                }else if(ff4[i0,j0]<=0){
                                    leftvolume=0;
                                }else{
                                    leftvolume=ff4[i0,j0]/2.0;
                                }
                            }
                                c_temp[i0+1][j0+1] = rightvolume + leftvolume;
                                // printf("leftvolume=%g,rightvolume=%g c=%g\n",leftvolume,rightvolume,rightvolume + leftvolume);  
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
          if((complete1.y[]==3 && modphase1.y[]>HUGE/2.0) || complete1.y[]==5){
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
                        double n_temp0_abs=0;
                      n_temp0_abs = sqrt(n_temp0.x*n_temp0.x+n_temp0.y*n_temp0.y);
                       if(fabs(n_temp0_abs)>0){  
                          b_left.y = 0.0;
                          leftvolume=rectangle_fraction(n_temp,alpha_temp,a_left,b_left)*(b_left.y-a_left.y);        
                       }else{
                          if(ff4[i0,j0-1]>=1){
                                    leftvolume=0.5;
                                }else if(ff4[i0,j0-1]<=0){
                                    leftvolume=0;
                                }else{
                                    leftvolume=ff4[i0,j0-1]/2.0;
                                }
                       }
                       double n_temp_abs=0;
                        n_temp_abs = sqrt(n_temp.x*n_temp.x+n_temp.y*n_temp.y);
                       if(fabs(n_temp_abs)>0){
                          a_right.y = 0.0;
                          rightvolume=rectangle_fraction(n_temp0,alpha_temp0,a_right,b_right)*(b_right.y-a_right.y);  
                       }else{
                          if(ff4[i0,j0]>=1){
                                    leftvolume=0.5;
                                }else if(ff4[i0,j0]<=0){
                                    leftvolume=0;
                                }else{
                                    leftvolume=ff4[i0,j0]/2.0;
                                }
                       }          
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
     boundary({modphase1,modphase0});
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



// // // // extern vector complete1_show;


// // // // // void get_modphase01_3(scalar ff4, face vector modphase1, face vector modphase0, face vector modphasexyl, face vector modphaser,vector nn4,scalar area4){
// // // //   void get_modphase01_3(scalar ff4, face vector modphase1, face vector modphase0, face vector modphasexyl, face vector modphaser,vector nn4){


// // // // //check mod//initial
// // // // // foreach(){
// // // // //     foreach_dimension(){
// // // // //       smallmodl.x[]=-10.0;
// // // // //       bigmodl.x[]=-10.0;
// // // // //       smallmodg.x[]=-10.0;
// // // // //       bigmodg.x[]=-10.0;

// // // // //     }
// // // // // }
// // // // ///////////////
    
    
// // // //     // scalar ff4=q.ff4;
// // // //     // face vector modphase1=q.modphase1;
// // // //     // face vector modphase0=q.modphase0;
// // // //     // scalar Tl=q.Tl;
// // // //     vector hhh[];
// // // //     heights(ff4,hhh);
// // // //     // vector nn4[];
// // // //     scalar alpha4[];
// // // //     scalar area4[];
// // // //      foreach_dimension()
// // // //         nn4.x.refine = nn4.x.prolongation = refine_injection;

// // // //     // vector hhh=q.hhh;
// // // //     // foreach(){
// // // //     //    ff4[]=ff[];
// // // //     // }   
// // // //     // ff4.restriction = restriction_volume_average;
// // // //     // ff4.restriction = restriction_volume_average;
// // // //     ff4.refine = ff4.prolongation = fraction_refine;
    
// // // //     boundary({hhh});
    
// // // //     restriction({ff4});

// // // // //taken from restruction
// // // // foreach() {
// // // //     area4[] = arealg[];
// // // //     if (ff4[] <= 0. || ff4[] >= 1.) {
// // // //       alpha4[] = 0.;
// // // //       // area4[]=0.0;
// // // //       foreach_dimension()
// // // // 	       nn4.x[] = 0.;
// // // //       }else{
// // // //       coord m = interface_normal3 (point, ff4, hhh);
// // // //       foreach_dimension()
// // // // 	        nn4.x[] = m.x;
// // // //       // alpha4[] = plane_alpha (ff4[], m);
// // // //       // coord temp;
// // // //       // double area = plane_area_center (m, alpha4[], &temp);
// // // //       // if (metric_embed_factor)  
// // // //       //       area *= metric_embed_factor (point, temp);
// // // //       // area4[] = area;
// // // //     }
// // // //   }

// // // // #if TREE
// // // //   foreach_dimension()
// // // //     nn4.x.refine = nn4.x.prolongation = refine_injection;

// // // //   alpha4.n = nn4;
// // // //   alpha4.refine = alpha4.prolongation = alpha_refine;
// // // // #endif
  

// // // //   //double lim_cut = 1e-12;
// // // //   double lim_cut = 1e-2;
// // // //   restriction({ff4,nn4,alpha4});
// // // //   face vector complete1[];
// // // //   // for(int l=0;l<=depth();l++){
// // // //   //      foreach_level(l){
// // // //   //         foreach_dimension(){
// // // //   //             complete1.x[] = 0;
// // // //   //             modphase1.x[] = HUGE;
// // // //   //             modphase0.x[] = HUGE;
// // // //   //             Point neib = neighborp(1);
// // // //   //             if(is_boundary(neighbor(1))){
// // // //   //                  foreach_neighbor(1){
// // // //   //                     if(point.i==neib.i && point.j==neib.j && point.k==neib.k){
// // // //   //                        modphase1.x[] = HUGE;
// // // //   //                        modphase0.x[] = HUGE;
// // // //   //                        complete1.x[] = 0;
// // // //   //                     }
// // // //   //                  }
// // // //   //             }
// // // //   //         }
// // // //   //      }
// // // //   // }

// // // //   foreach_face(){
// // // //           complete1.x[] = 0;
// // // //           modphase1.x[] = HUGE;
// // // //           modphase0.x[] = HUGE;
// // // //           modphasexyl.x[] = HUGE;
// // // //           modphasexyr.x[] = HUGE;
// // // //   }
// // // //   // for(int l=0;l<=depth();l++){
// // // //   //    foreach_level(l){
// // // //   //        foreach_dimension(){
// // // //           foreach_face(){
// // // //                if((ff4[-1]==0.0 && ff4[]==1.0) || (ff4[-1]==1.0 && ff4[]==0.0)){
// // // //                      modphase1.x[]=0.5;
// // // //                      modphase0.x[]=0.5;
// // // //                      modphasexyl.x[]=0.5;
// // // //                      modphasexyr.x[]=0.5;
// // // //                      complete1.x[] = 6;
// // // //                 }else if( (is_phase(ff4[-1]) && !is_phase(ff4[])) || (!is_phase(ff4[-1]) && is_phase(ff4[]))){ // 
// // // //                     complete1.x[]=3;
// // // // 		                 if((is_phase(ff4[-1]) && !is_phase(ff4[])) && orientation(hhh.x[])==0){ // get from height,(1 0)           
// // // //                         if(height(hhh.x[])>-1.0 && height(hhh.x[])<0.0){
// // // // 		                      modphase0.x[]=fabs(height(hhh.x[]));
// // // //                           modphase1.x[]=1.0-modphase0.x[];
// // // //                           modphasexyl.x[] = modphase1.x[];
// // // //                           modphasexyr.x[] = 1 -modphasexyl.x[];
// // // //                           complete1.x[]=2;
// // // //                          }
                       
// // // // 		                 }else if((!is_phase(ff4[-1]) && is_phase(ff4[])) && orientation(hhh.x[])==1){ // 0 1
// // // // 	                       if(height(hhh.x[])<0.0 && height(hhh.x[]>-1.0)){
// // // //                            modphase1.x[]=fabs(height(hhh.x[]));
// // // //                            modphase0.x[]=1.0-modphase1.x[];
// // // //                            modphasexyl.x[] = modphase0.x[];
// // // //                            modphasexyr.x[] = 1 - modphasexyl.x[];
// // // //                            complete1.x[]=2;
// // // //                           }
// // // // 		                 }
// // // //                  if(complete1.x[]==3){ // get from middle cell
// // // // 		                 if((is_phase(ff4[-1]) && !is_phase(ff4[])) && orientation(hhh.x[-1])==0){ // get from height,(1 0)           
// // // //                         if(height(hhh.x[-1])>0.0 && height(hhh.x[-1])<1.0){
// // // // 		                      modphase1.x[]=fabs(height(hhh.x[-1]));
// // // //                           modphase0.x[]=1.0-modphase1.x[];

// // // //                           modphasexyl.x[] = modphase1.x[];
// // // //                           modphasexyr.x[] = 1 -modphasexyl.x[];
// // // //                           complete1.x[]=2;
// // // //                          }
                       
// // // // 		                 }else if((!is_phase(ff4[-1]) && is_phase(ff4[])) && orientation(hhh.x[-1])==1){ // 0 1
// // // // 	                       if(height(hhh.x[-1])>0.0 && height(hhh.x[-1]<1.0)){
// // // //                            modphase0.x[]=fabs(height(hhh.x[-1]));
// // // //                            modphase1.x[]=1.0-modphase0.x[];

// // // //                            modphasexyl.x[] = modphase0.x[];
// // // //                            modphasexyr.x[] = 1 - modphasexyl.x[];
// // // //                            complete1.x[]=2;
// // // //                           }
// // // // 		                 }
// // // //                    }
// // // //                   // if(complete1.x[]==2 && (modphase1.x[]<lim_cut)){
// // // //                   //            modphase1.x[] = lim_cut;
// // // //                   //            modphase0.x[] = 1.0 - modphase1.x[];
// // // //                   //            complete1.x[] = 5;
// // // //                   //   }else if(complete1.x[]==2 && (modphase1.x[]>1.0-lim_cut)){
// // // //                   //            modphase1.x[] = 1.0 - lim_cut;
// // // //                   //            modphase0.x[] = 1.0 - modphase1.x[];
// // // //                   //            complete1.x[] = 5;
// // // //                   //  } 
                  
// // // //                   //  if(complete1.x[]==2 && (modphasexyl.x[]<lim_cut)){
// // // //                   //            modphase1.x[] = lim_cut;
// // // //                   //            modphase0.x[] = 1.0 - modphase1.x[];
// // // //                   //            complete1.x[] = 5;
// // // //                   //   }else if(complete1.x[]==2 && (modphasexyl.x[]>1.0-lim_cut)){
// // // //                   //            modphase1.x[] = 1.0 - lim_cut;
// // // //                   //            modphase0.x[] = 1.0 - modphase1.x[];
// // // //                   //            complete1.x[] = 5;
// // // //                   //  } 
// // // //                   //  modphase1.x[] = lim_cut;
// // // //                             //  modphase0.x[] = 1.0 - modphase1.x[];
// // // //                             if(modphase1.x[]<=0 || modphase1.x[]>=1){
// // // //                                 modphase1.x[]=1;
// // // //                                 modphase0.x[]=1;
// // // //                                 modphasexyl.x[]=1;
// // // //                                 modphasexyr.x[]=1;
// // // //                                 complete1.x[] = 5;
// // // //                             }
// // // // 	           }else{
// // // //                   modphase1.x[] = 1.0;
// // // //                   modphase0.x[] = 1.0;

// // // //                   modphasexyl.x[] = 1;
// // // //                   modphasexyr.x[] = 1;
// // // //                   complete1.x[] = 7;
// // // // 	          }
             
                            
                             
                    
// // // //         }
// // // //     foreach_face (x){
// // // //           if((complete1.x[]==3 && modphase1.x[]>HUGE/2.0) || complete1.x[]==5){
// // // //                double c_temp[3][3];
// // // //                for(int i0=-1;i0<2;i0++){
// // // //                   for(int j0=-1;j0<2;j0++){
// // // //                             //coord leftvolume,rightvolume;  //different from linear2.h, rifhtvolume is the right volume of -1; left volume if left volume of 0;
// // // //                             double leftvolume,rightvolume;
// // // //                             coord a_left,b_left;
// // // //                             coord a_right,b_right;
// // // //                             coord n_temp;
// // // //                             n_temp.x = nn4.x[i0,j0], n_temp.y = nn4.y[i0,j0];

// // // //                             coord n_temp0;
// // // //                             n_temp0.x = nn4.x[i0-1,j0], n_temp0.y = nn4.y[i0-1,j0];
// // // //                             double alpha_temp = alpha4[i0,j0];
// // // //                             double alpha_temp0 = alpha4[i0-1,j0];
// // // //                             //foreach_dimension(){
                                
// // // //                                 a_left=(coord){-0.5,-0.5};
// // // //                                 b_left=(coord){0.5,0.5};
                                
                                    
// // // //                                 a_right=(coord){-0.5,-0.5};
// // // //                                 b_right=(coord){0.5,0.5};   
// // // //                                double n_temp0_abs=0;
// // // //                             n_temp0_abs = sqrt(n_temp0.x*n_temp0.x+n_temp0.y*n_temp0.y);
// // // //                             if(fabs(n_temp0_abs)>0){
// // // //                                 b_left.x = 0.0;
// // // //                                 leftvolume=rectangle_fraction(n_temp,alpha_temp,a_left,b_left)*(b_left.x-a_left.x);        
// // // //                             }else{
// // // //                                 if(ff4[i0-1,j0]>=1){
// // // //                                     leftvolume=0.5;
// // // //                                 }else if(ff4[i0-1,j0]<=0){
// // // //                                     leftvolume=0;
// // // //                                 }else{
// // // //                                     leftvolume=ff4[i0-1,j0]/2.0;
// // // //                                 }
// // // //                             }
// // // //                             double n_temp_abs=0;
// // // //                             n_temp_abs = sqrt(n_temp.x*n_temp.x+n_temp.y*n_temp.y);
// // // //                             if(fabs(n_temp_abs)>0){
// // // //                                 a_right.x = 0.0;
// // // //                                 rightvolume=rectangle_fraction(n_temp0,alpha_temp0,a_right,b_right)*(b_right.x-a_right.x);  
// // // //                             }else{
// // // //                                if(ff4[i0,j0]>=1){
// // // //                                     leftvolume=0.5;
// // // //                                 }else if(ff4[i0,j0]<=0){
// // // //                                     leftvolume=0;
// // // //                                 }else{
// // // //                                     leftvolume=ff4[i0,j0]/2.0;
// // // //                                 }
// // // //                             }
                                                  
// // // //                                 c_temp[i0+1][j0+1] = rightvolume + leftvolume;
// // // //                           //}// foreach_dimension
// // // //                   } //for j0
// // // //                 } //for i0
// // // //                 //I am here
// // // //                 coord n_middle = mycs2(c_temp);//in myc2d.h //mycs3(c_temp);  //get normal from middle volume distribution
// // // //                 double c_middle = c_temp[0+1][0+1];
// // // //                 double alpha_middle = plane_alpha(c_middle,n_middle);
// // // //                 double test=0;
                
// // // //                 if(fabs(n_middle.x) > 1e-12){ //n.x x + n.y y = alpha y=0 kan x
// // // // 	           test =  alpha_middle/n_middle.x - (-0.5);
// // // // 		   // fprintf(fp18,"nnx=%g nny=%g test=%g\n",nn.x,nn.y,test);
// // // //                 }
// // // //                 if(c_middle<=0 || c_middle>=1){
// // // //                     modphasexyl.x[] = 1;
// // // //                     modphasexyr.x[] = 1;

// // // //                      modphasexyl.z[] = 1;
// // // //                     modphasexyr.x[] = 1;

// // // //                     complete1.x[] = 8;     
// // // //                 }else{
// // // //                     if(test<=0 || test>=1){
// // // //                         modphase1.x[]=1;
// // // //                         modphase0.x[]=1;  

// // // //                         modphasexyl.x[] = 1;
// // // //                         modphasexyr.x[] = 1;
// // // //                         complete1.x[] = 8;     
// // // //                     }else{
// // // //                         modphase1.x[]=test;
// // // //                         modphase0.x[]=1.0-modphase1.x[];  

// // // //                         if(n_middle.x>0){
// // // //                           modphasexyl.x[] = test;
// // // //                           modphasexyr.x[] = 1-test;
// // // //                           complete1.x[] = 4;     
// // // //                         }else if(n_middle.x<0){
// // // //                           modphasexyl.x[] = 1-test;
// // // //                           modphasexyr.x[] = test;
// // // //                           complete1.x[] = 4;     
// // // //                         }else{ //==0
// // // //                           modphase1.x[]=1;
// // // //                           modphase0.x[]=1;  

// // // //                           modphasexyl.x[] = 1;
// // // //                           modphasexyr.x[] = 1;
// // // //                           complete1.x[] = 8;     
// // // //                         }
// // // //                     }
// // // //                 }



// // // //            //     double lim_cut = 1e-13;
// // // //                 // if(test>lim_cut && test<1-lim_cut)
// // // //                 //   {
// // // //                 //         test = test;
// // // //                 //   }else if(test>1.0-lim_cut){
// // // //                 //         test = 1.0 - lim_cut;
// // // //                 //   }else if(test<lim_cut){
// // // //                 //         test = lim_cut;
// // // //                 //   }

// // // //                 // modphase1.x[]=test;
// // // //                 // modphase0.x[]=1.0-modphase1.x[];  
// // // //                 //complete.x = 4; 
                      
// // // //    //
// // // //           }//complete.x==3
// // // //           if(complete1.x[]==8){
// // // //               coord n_temp;
// // // //               n_temp.x = nn4.x[], n_temp.y = nn4.y[];

// // // //               coord n_temp0;
// // // //               n_temp0.x = nn4.x[-1,0], n_temp0.y = nn4.y[-1,0];
// // // //               double alpha_temp = alpha4[];
// // // //               double alpha_temp0 = alpha4[-1,0];
              
// // // //                double test0=0;
                
// // // //                 if(fabs(n_temp0.x) > 1e-12){ //n.x x + n.y y = alpha y=0 kan x
// // // // 	                 test0 =  alpha_temp0/n_temp0.x;

// // // //                         if(n_temp0.x>0){
// // // //                           if(test0>=0.5){
// // // //                               test0=0.5; //left
// // // //                               test0=0.5-test0; //right
// // // //                           }else if(test0<=0){
// // // //                               test0=0; //left
// // // //                               test0=0.5-0; //right
// // // //                           }else{
// // // //                               test0=test0;//left
// // // //                               test0=0.5-test0;//right
// // // //                           }
// // // //                           // modphasexyl.x[] = test;
// // // //                           // modphasexyr.x[] = 1-test;
// // // //                           // complete1.x[] = 4;     
// // // //                         }else if(n_temp0.x<0){  
// // // //                            if(test0>=0.5){
// // // //                               test0=0.5; //right
// // // //                           }else if(test0<=0){
// // // //                               test0=0; //right
// // // //                           }else{
// // // //                               test0=0.5-test0;//left
// // // //                               test0= 0.5-test0;//right
// // // //                           }
// // // //                           // modphasexyl.x[] = 1-test;
// // // //                           // modphasexyr.x[] = test;
// // // //                           // complete1.x[] = 4  
// // // //                         }else{ //==0
// // // //                           if(ff4[-1,0]<=0){
// // // //                               // left 0->1
// // // //                               test0=-1;//() gas;
// // // //                           }else if(ff4[-1,0]>=1){
// // // //                               test0=-2;// liquid
// // // //                           }    
// // // //                         }
                    
// // // //                 }

// // // //             double test=0;
                
// // // //                 if(fabs(n_temp.x) > 1e-12){ //n.x x + n.y y = alpha y=0 kan x
// // // // 	                 test =  alpha_temp/n_temp.x;

// // // //                        if(n_temp.x>0){
// // // //                           if(test>=0){
// // // //                               test=0.5; //left
// // // //                               // test0=0.5-test0; //right
// // // //                           }else if(test<=-0.5){
// // // //                               test=0; //left
// // // //                               // test0=0.5-0; //right
// // // //                           }else{
// // // //                               test=test-(-0.5);//left
// // // //                               // test0=0.5-test0;//right
// // // //                           }
// // // //                           // modphasexyl.x[] = test;
// // // //                           // modphasexyr.x[] = 1-test;
// // // //                           // complete1.x[] = 4;     
// // // //                         }else if(n_temp.x<0){
// // // //                           if(test>=0.5){
// // // //                               test=0; //left
// // // //                           }else if(test<=0){
// // // //                               test=0.5; //left
// // // //                           }else{
// // // //                               test=0.5-test;//left
// // // //                               // test0= 0.5-test0//right
// // // //                           }
// // // //                           // modphasexyl.x[] = 1-test;
// // // //                           // modphasexyr.x[] = test;
// // // //                           // complete1.x[] = 4;         
// // // //                         }else{ //==0
// // // //                           if(ff4[]<=0){
// // // //                               // left 0->1
// // // //                               test=-1;//() gas;
// // // //                           }else if(ff4[]>=1){
// // // //                               test=-2;// liquid
// // // //                           }    
// // // //                         }
                    
// // // //                 }
// // // //                 if(test0==(-1) && test==(-1)){
// // // //                     modphasexyl.x[]=1;
// // // //                     modphasexyr.x[]=1;
// // // //                 }else if(test0==(-2) && test==(-2)){
// // // //                     modphasexyl.x[]=1;
// // // //                     modphasexyr.x[]=1;
// // // //                 }else if((test0==(-2) && test==(-1))||(test0==(-1) && test==(-2))){
// // // //                     modphasexyl.x[]=0.5;
// // // //                     modphasexyr.x[]=0.5;
// // // //                 }else {
// // // //                     if(test0==-1 || test0==-2){
// // // //                         test0=0.5;
// // // //                     }
// // // //                     if(test==-1 || test==-2){
// // // //                         test=0.5;
// // // //                     }
// // // //                     modphasexyl.x[]=test0+test;
// // // //                     if(modphasexyl.x[]<=0 || modphasexyl.x[]>=1){
// // // //                         modphasexyl.x[]=1;
// // // //                         modphasexyr.x[]=1;
// // // //                     }else{
// // // //                         modphasexyr.x[]=1.0-modphasexyl.x[];
// // // //                     }
// // // //                 }
// // // //                 complete1.x[]=9;
// // // //           }
// // // //   } //foreach_face
// // // //     foreach_face (y){
// // // //           if((complete1.y[]==3 && modphase1.y[]>HUGE/2.0) || (complete1.y[]==5) ){
// // // //                double c_temp[3][3];
// // // //                for(int i0=-1;i0<2;i0++){
// // // //                   for(int j0=-1;j0<2;j0++){
// // // //                           //coord leftvolume,rightvolume;  //different from linear2.h, rifhtvolume is the right volume of -1; left volume if left volume of 0;
// // // //                           double leftvolume,rightvolume;
// // // //                           coord a_left,b_left;
// // // //                           coord a_right,b_right;
// // // //                               //coord n_temp= interface_normal (point, ff);
// // // //                           coord n_temp;
// // // //                           n_temp.x = nn4.x[i0,j0], n_temp.y = nn4.y[i0,j0];
// // // //                           coord n_temp0;
// // // //                           n_temp0.x = nn4.x[i0,j0-1], n_temp0.y = nn4.y[i0,j0-1];
// // // //                               //double alpha_temp = plane_alpha (ff[], n_temp);
// // // //                           double alpha_temp = alpha4[i0,j0];
// // // //                           double alpha_temp0 = alpha4[i0,j0-1];
// // // //                               //foreach_dimension(){
                              
// // // //                           a_left=(coord){-0.5,-0.5};
// // // //                           b_left=(coord){0.5,0.5};
                          
                              
// // // //                           a_right=(coord){-0.5,-0.5};
// // // //                           b_right=(coord){0.5,0.5};
// // // //                         double n_temp0_abs=0;
// // // //                       n_temp0_abs = sqrt(n_temp0.x*n_temp0.x+n_temp0.y*n_temp0.y);
// // // //                        if(fabs(n_temp0_abs)>0){  
// // // //                           b_left.y = 0.0;
// // // //                           leftvolume=rectangle_fraction(n_temp,alpha_temp,a_left,b_left)*(b_left.y-a_left.y);        
// // // //                        }else{
// // // //                           if(ff4[i0,j0-1]>=1){
// // // //                                     leftvolume=0.5;
// // // //                                 }else if(ff4[i0,j0-1]<=0){
// // // //                                     leftvolume=0;
// // // //                                 }else{
// // // //                                     leftvolume=ff4[i0,j0-1]/2.0;
// // // //                                 }
// // // //                        }
// // // //                        double n_temp_abs=0;
// // // //                         n_temp_abs = sqrt(n_temp.x*n_temp.x+n_temp.y*n_temp.y);
// // // //                        if(fabs(n_temp_abs)>0){
// // // //                           a_right.y = 0.0;
// // // //                           rightvolume=rectangle_fraction(n_temp0,alpha_temp0,a_right,b_right)*(b_right.y-a_right.y);  
// // // //                        }else{
// // // //                           if(ff4[i0,j0]>=1){
// // // //                                     leftvolume=0.5;
// // // //                                 }else if(ff4[i0,j0]<=0){
// // // //                                     leftvolume=0;
// // // //                                 }else{
// // // //                                     leftvolume=ff4[i0,j0]/2.0;
// // // //                                 }
// // // //                        }                   
// // // //                           c_temp[i0+1][j0+1] = rightvolume + leftvolume;
// // // //                             //}// foreach_dimension
// // // //                   } //for j0
// // // //                 } //for i0
// // // //                 coord n_middle = mycs2(c_temp);//in myc2d.h //mycs3(c_temp);  //get normal from middle volume distribution
// // // //                 double c_middle = c_temp[0+1][0+1];
// // // //                 double alpha_middle = plane_alpha(c_middle,n_middle);
// // // //                 double test=0;;
// // // //                 if(fabs(n_middle.y) > 1e-12){
// // // // 	                   test =  alpha_middle/n_middle.y - (-0.5);
// // // // 		                 // fprintf(fp18,"nnx=%g nny=%g test=%g\n",nn.x,nn.y,test);
// // // //                 } 
// // // //              //   double lim_cut = 1e-13;

// // // //                 if(c_middle<=0 || c_middle>=1){
// // // //                     modphase1.y[]=1;
// // // //                     modphase0.y[]=1;

// // // //                     modphasexyl.y[] = 1;
// // // //                     modphasexyr.y[] = 1;

// // // //                     complete1.y[]=8;
// // // //                 }else{
// // // //                     if(test<=0 || test>=1){
// // // //                         modphase1.y[]=1;
// // // //                         modphase0.y[]=1;

// // // //                         modphasexyl.y[] = 1;
// // // //                         modphasexyr.y[] = 1;

// // // //                         complete1.y[]=8;
// // // //                     }else{
// // // //                        modphase1.y[]=test;
// // // //                        modphase0.y[]=1.0-modphase1.y[];

// // // //                        if(n_middle.y>0){
// // // //                           modphasexyl.y[] = test;
// // // //                           modphasexyr.y[] = 1-test;

// // // //                            complete1.y[] = 4;   
// // // //                         }else if(n_middle.y<0){
// // // //                           modphasexyl.y[] = 1-test;
// // // //                           modphasexyr.y[] = test;
// // // //                            complete1.y[] = 4;   
// // // //                         }else{ //==0
// // // //                           modphase1.y[]=1;
// // // //                           modphase0.y[]=1;


// // // //                           modphasexyl.y[] = 1;
// // // //                           modphasexyr.y[] = 1;
// // // //                            complete1.y[] = 8;   
// // // //                         }
// // // //                     }
// // // //                 }



// // // //                 // if(test>lim_cut && test<1-lim_cut)
// // // //                 //   {
// // // //                 //         test = test;
// // // //                 //   }else if(test>1.0-lim_cut){
// // // //                 //         test = 1.0 - lim_cut;
// // // //                 //   }else if(test<lim_cut){
// // // //                 //         test = lim_cut;
// // // //                 //   }

// // // //                 // modphase1.y[]=test;
// // // //                 // modphase0.y[]=1.0-modphase1.y[];
// // // //                // complete.y = 4; 
                           
// // // //    //           
// // // //           }// complete.y==3
// // // //          if(complete1.y[]==8){
// // // //               coord n_temp;
// // // //               n_temp.x = nn4.x[], n_temp.y = nn4.y[];

// // // //               coord n_temp0;
// // // //               n_temp0.x = nn4.x[0,-1], n_temp0.y = nn4.y[0,-1];
// // // //               double alpha_temp = alpha4[];
// // // //               double alpha_temp0 = alpha4[0,-1];
              
// // // //                double test0=0;
                
// // // //                 if(fabs(n_temp0.y) > 1e-12){ //n.x x + n.y y = alpha y=0 kan x
// // // // 	                 test0 =  alpha_temp0/n_temp0.y;

// // // //                         if(n_temp0.y>0){
// // // //                           if(test0>=0.5){
// // // //                               test0=0.5; //left
// // // //                               test0=0.5-test0; //right
// // // //                           }else if(test0<=0){
// // // //                               test0=0; //left
// // // //                               test0=0.5-0; //right
// // // //                           }else{
// // // //                               test0=test0;//left
// // // //                               test0=0.5-test0;//right
// // // //                           }
// // // //                           // modphasexyl.x[] = test;
// // // //                           // modphasexyr.x[] = 1-test;
// // // //                           // complete1.x[] = 4;     
// // // //                         }else if(n_temp0.y<0){
// // // //                            if(test0>=0.5){
// // // //                               test0=0.5; //right
// // // //                           }else if(test0<=0){
// // // //                               test0=0; //right
// // // //                           }else{
// // // //                               test0=0.5-test0;//left
// // // //                               test0= 0.5-test0;//right
// // // //                           }
// // // //                           // modphasexyl.x[] = 1-test;
// // // //                           // modphasexyr.x[] = test;
// // // //                           // complete1.x[] = 4;     
// // // //                         }else{ //==0
// // // //                           if(ff4[0,-1]<=0){
// // // //                               // left 0->1
// // // //                               test0=-1;//() gas;
// // // //                           }else if(ff4[0,-1]>=1){
// // // //                               test0=-2;// liquid
// // // //                           }    
// // // //                         }
                    
// // // //                 }

// // // //             double test=0;
                
// // // //                 if(fabs(n_temp.y) > 1e-12){ //n.x x + n.y y = alpha y=0 kan x
// // // // 	                 test =  alpha_temp/n_temp.y;

// // // //                        if(n_temp.y>0){
// // // //                           if(test>=0){
// // // //                               test=0.5; //left
// // // //                               // test0=0.5-test0; //right
// // // //                           }else if(test<=-0.5){
// // // //                               test=0; //left
// // // //                               // test0=0.5-0; //right
// // // //                           }else{
// // // //                               test=test-(-0.5);//left
// // // //                               // test0=0.5-test0;//right
// // // //                           }
// // // //                           // modphasexyl.x[] = test;
// // // //                           // modphasexyr.x[] = 1-test;
// // // //                           // complete1.x[] = 4;     
// // // //                         }else if(n_temp.y<0){
// // // //                           if(test>=0.5){
// // // //                               test=0; //left
// // // //                           }else if(test<=0){
// // // //                               test=0.5; //left
// // // //                           }else{
// // // //                               test=0.5-test;//left
// // // //                               // test0= 0.5-test0//right
// // // //                           }
// // // //                           // modphasexyl.x[] = 1-test;
// // // //                           // modphasexyr.x[] = test;
// // // //                           // complete1.x[] = 4;     
// // // //                         }else{ //==0
// // // //                           if(ff4[]<=0){
// // // //                               // left 0->1
// // // //                               test=-1;//() gas;
// // // //                           }else if(ff4[]>=1){
// // // //                               test=-2;// liquid
// // // //                           }    
// // // //                         }
                    
// // // //                 }
// // // //                 if(test0==(-1) && test==(-1)){
// // // //                     modphasexyl.y[]=1;
// // // //                     modphasexyr.y[]=1;
// // // //                 }else if(test0==(-2) && test==(-2)){
// // // //                     modphasexyl.y[]=1;
// // // //                     modphasexyr.y[]=1;
// // // //                 }else if((test0==(-2) && test==(-1))||(test0==(-1) && test==(-2))){
// // // //                     modphasexyl.y[]=0.5;
// // // //                     modphasexyr.y[]=0.5;
// // // //                 }else {
// // // //                     if(test0==-1 || test0==-2){
// // // //                         test0=0.5;
// // // //                     }
// // // //                     if(test==-1 || test==-2){
// // // //                         test=0.5;
// // // //                     }
// // // //                     modphasexyl.y[]=test0+test;
// // // //                     if(modphasexyl.y[]<=0 || modphasexyl.y[]>=1){
// // // //                         modphasexyl.y[]=1;
// // // //                         modphasexyr.y[]=1;
// // // //                     }else{
// // // //                         modphasexyr.y[]=1.0-modphasexyl.y[];
// // // //                     }
// // // //                 }
// // // //                 complete1.y[]=9;
// // // //           }
          
// // // //   } 
  
// // // //   //foreach_face (y)       
// // // //         //  foreach_dimension(){
// // // //         //     if(is_boundary(neighbor(1))){
// // // //         //         modphase1.x[]=
// // // //         //     }
// // // //         //  }


// // // //    //set boundary for mod
// // // //     foreach_face() {
         
// // // //     }
// // // //      	  foreach_face(){
// // // //                      if(fabs(modphase0.x[])<lim_cut){
// // // //                           modphase0.x[]=lim_cut;
// // // //                           modphase1.x[]=1.0-modphase0.x[];
// // // //                      }else if(fabs(modphase1.x[])<lim_cut){
// // // //                           modphase1.x[]=lim_cut;
// // // //                           modphase0.x[]=1.0-modphase1.x[];
// // // //                      }
// // // //                      if(modphase0.x[]>1){
// // // //                           modphase0.x[] = 1.0;
// // // //                      }
// // // //                      if(modphase1.x[]>1){
// // // //                           modphase1.x[] = 1.0;
// // // //                      }
// // // //          }

// // // //     foreach_dimension(){
// // // //       modphase0.x.restriction = face_rejection2_x ;
// // // //       //   modphase0.x.prolongation = mod_prol_x;
// // // //          modphase1.x.restriction = face_rejection2_x ;
// // // //           modphasexyl.x.restriction = face_rejection2_x ;
// // // //           modphasexyr.x.restriction = face_rejection2_x ;
// // // //         // modphase0.x.prolongation = mod_prol_x;
// // // //     }
// // // //        restriction({modphase0.x, modphase0.y, modphase1.x,modphase1.y});
// // // //        restriction({modphasexyl.x, modphasexyl.y, modphasexyr.x,modphasexyr.y});
// // // //      boundary({modphase1,modphase0});
// // // //      boundary({modphasexyl,modphasexyr});
// // // //     //  char name101[80];
// // // // 	  //  sprintf(name101,"mass_record5-pid%d.dat",pid());
// // // // 	  //   FILE * fp101 = fopen(name101,"w");
// // // //     //   foreach_face(){
// // // //     //      // if(fabs(modphase1.x[])<=1.0){
// // // //     //             fprintf(fp101,"%g %g %g %g %g\n",x,y,z,modphase1.x[],modphase0.x[]);
// // // //     //       //}
// // // //     //   }
// // // //     // fclose(fp101);
// // // //     //   MPI_Barrier(MPI_COMM_WORLD);
// // // //     //   if(pid()==0){
// // // //     //       char command1[150];
// // // //     //       sprintf(command1, "LC_ALL=C cat mass_record5-pid*.dat > outfacets/mass_record5-%g",t);
// // // //     //       system(command1);

// // // //     //       char command7[150];
// // // //     //       sprintf(command7, "LC_ALL=C rm -rf mass_record5-pid*.dat");
// // // //     //       system(command7);
// // // //     //   }

// // // //     foreach(){
// // // //       complete1_show.x[] = complete1.x[];
// // // //       complete1_show.y[] = complete1.y[];
// // // //     }
  
// // // //   }
///for masstr[] source_pc[] 
//scalar masstr[];
extern scalar masstr;
extern scalar source_pc;
extern scalar source_pc2;
extern scalar vtr;
extern double Tkg,Tkl,Trhog,Trhol,Tcpg,Tcpl,hfg;
extern double source_total,total_area;

extern scalar phase0Tgrad,phase1Tgrad;
extern scalar aiml,aimg;
//extern scalar source_ff1,source_ff2,ff_old1;
//#if TGRAD_LEON
#if 1
     #include "Tgrad-leon.h"
#endif 

#include "./my-tag.h"

void mass_transfer_rate(){
   //scalar masstr[];
   
   source_total= 0.0;
   total_area=0.0;
  // scalar phase0Tgrad[],phase1Tgrad[];
//#if TGRAD_LEON
#if 1
     Tgrad_leon(ff,phase0Tgrad,phase1Tgrad,aiml,aimg);
      // Tgrad_leon(ff,phase0Tgrad,phase1Tgrad);
      // printf("Trad_leon\n");
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


  

  scalar d_tag[];
  double threshold = 1e-10;
  foreach()
    d_tag[] =  (topo_mask_s[]==0 && topo_mask[]==0)? (ff[] > threshold):0;
  int n_tag = tag (d_tag), size[n_tag];
   double total_source[n_tag];
 double area_source[n_tag];
  for (int i = 0; i < n_tag; i++){
    size[i] = 0;
    total_source[i]=0.0;
    area_source[i] = 0.0;
  }
  foreach (serial)
    if (d_tag[] > 0)
      size[((int) d_tag[]) - 1]++;
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, size, n_tag, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  // int minsize = pow (p.minsize ? p.minsize : 3, dimension);
  // foreach()
  //   if (d[] > 0 && size[((int) d[]) - 1] < minsize)
  //     f[] = p.bubbles;




   foreach(reduction(+:source_total) reduction(+:total_area) reduction(+:total_source[:n_tag]) 
   reduction(+:area_source[:n_tag])){
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
                  area = area*max(1e-20,(y+pp.y*Delta));
                //  source_pc2[] = 0.0; //source_pc[]*(y+pp.y*Delta); //20230603
                // printf("source_pc9999\n");
        #endif
        source_total = source_total + source_pc2[];
        if(d_tag[]>0){
          total_source[((int)d_tag[])-1]+= source_pc2[];
        }

        if(cs[]>0.0){
          if(d_tag[]>0){
            total_area = total_area + area;
            area_source[((int)d_tag[])-1]+=area;
          }
        }

       // source_ff1[] = 1; //cells contain sources;
   }

  //  if(t<0.000050){
    if(flag_average_source){
        foreach(){
            if(cs[]>0.0){
                if(d_tag[]>0){
                    double area=0.0;
                    coord pp = {0.,0.};
                    if (ff[] > EPS && ff[] < 1. - EPS) {
                        coord n = interface_normal (point, ff);
                          
                        double alpha = plane_alpha (ff[], n);
                              area = line_length_center(n, alpha, &pp);//plane_area_center (n, alpha, &p);
                        //area += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
                      }

                      #if AXI
                                area = area*max(1e-20,(y+pp.y*Delta));
                      #endif
                      // source_pc2[] = source_total*area/total_area;
                      
                        source_pc2[] = total_source[((int)d_tag[])-1]*area/area_source[((int)d_tag[])-1];
                        // source_pc[] = masstr[]*(1.0/Trhog - 1.0/Trhol)*area/Delta;
                        masstr[]=source_pc2[]/((1.0/Trhog - 1.0/Trhol)*area/Delta);
                  }
            }else{
              source_pc2[] = 0.0;
            }
        }
    }

  // }

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


void  mov_interface_dc(bool flag_topos_uf,bool flag_cant_smaller_than_half,double base, double over_half){

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
        if(!flag_topos_uf){
            ff[] = ff[] + deltac[];
        }else{
            if(!(topo_mask_s[]==0 && intersect_true[]==0)){
            // if(!(topo_mask_s[]==0)){
              ff[] = ff[] + deltac[];
            }
        }
        
	    
   }

if(flag_cant_smaller_than_half){ //distribute over deltac to tripple point
  double total_dc_add=0.0;
  int total_dc_number=0;
  int number_interface=0;
  double y_interface_triple=0.0;
  foreach(reduction(+:total_dc_add) reduction(+:number_interface) 
   reduction(+:y_interface_triple) reduction(+:total_dc_number)){
    if(topo_mask_s[]==0){
      bool flag=false;
              // foreach_neighbor(){
      foreach_neighbor(2){
          if(intersect_true[]==1){
                      flag=true;
          }
      }
        if(f_old[]<1.0 && f_old[]>0.0 && (!flag)){
            // if((ff[]<base+overhalf) && f_old[]>base+overhalf){
            if((ff[]<base+over_half) && deltac[]<0.0){
                  total_dc_add += deltac[]*y; // true volume = y*area*distant = y*deltac
                  total_dc_number +=1;
                  ff[] = f_old[];
            } 
        }
    }
    if(intersect_true[]==1){
        number_interface+=1;
        y_interface_triple+=y;
    }
  }

  if(total_dc_number>=1 && number_interface>=1){
    foreach(){
        if(intersect_true[]==1){ //total_dc_add distributed by area, and then translate to no y
            // ff[] = f_old[] + (total_dc_add*y/y_interface_triple)/y;
            ff[] = f_old[] + y/y_interface_triple*total_dc_add/y;
        }
    }
  }
}

for(scalar s in {T}){
      foreach(){
           double val_tot=0.0;
           double wei=0.0;
           double val1=0.0;
        if(energy_advecting_flag){
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
      }else{
        if(css_test3_n[]>0.5){
                    if(ff[] >= 0.5 && f_old[] < 0.5 ){ // 20230605
                        val_tot = Tsat00*ff[];
                        wei = wei + ff[];
                       
                    }else if(1.0-ff[] > 0.5 && 1.0-f_old[] <= 0.5){
                        val_tot = val_tot + Tsat00*(1.0-ff[]);
                        wei = wei + (1.0-ff[]);
                    }else{
                        val_tot = s[];
                        wei = 1.0;
                    }
                    val1 = val_tot/wei;
           }
          //  if(css_test3_n[]<1.0){ //css_test3[]>0.0
          //       s[] = val1*css_test3_n[] + T_solid[]*(1.0-css_test3_n[]);
          //   }else{
          //       s[] = val1;
          //   }
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
      //  if(facemask.x[]==1){
          if (level==level_interface){
              int correct1 =0;

//  correct = ( (topo_mask(i,j,k)==-2 .and. topo_mask(i+i0,j+j0,k+k0)==-1) &
//          .or. (topo_mask(i,j,k)==-1) &
//          .or. (topo_mask(i,j,k)==0 .and. &
//          (topo_mask(i+i0,j+j0,k+k0)==0 .or. topo_mask(i+i0,j+j0,k+k0)==-1) ) )

              // if((topo_mask[-1]<1 && topo_mask[]<1) && (topo_mask[]==0 || topo_mask[-1]==0)  ){
                // //  if((topo_mask[-1]==-2 && topo_mask[]==-1) || 
                // //  (topo_mask[-1]==-1) || 
                // //  (topo_mask[-1]==0 && (topo_mask[-1]==0 || topo_mask[]==-1))  ){
                    if(flag_topos_advect_uf){
                        if(topo_mask_s[]==0 || topo_mask_s[-1]==0){ 
                            ulf.x[] = uf.x[];
                        }else{
                            ulf.x[] = uf.x[] + usf.x[];
                        }
                        
                    }else{
                        ulf.x[] = uf.x[] + usf.x[];
                    }
                      
                // // }else{
                // //       ulf.x[] = uf.x[];
                // // }
          }else{
              // ulf.x[] = uf.x[];
               ulf.x[] = uf.x[] + usf.x[];
          }

        // }
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

static void embed_fraction_refine_css_test (Point point, scalar css_test)
//void embed_fraction_refine_css_test (Point point, scalar css_test)
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

    // coord n = facet_normal (point, css_test, fss_test);
    coord n = mycs (point, ff);

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
static void embed_face_fraction_fss_test_refine_x (Point point, scalar s)
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

     // coord n = facet_normal (point, css_test, fss_test);
    coord n = mycs (point, ff);
    // double alpha = plane_alpha (css_test[], n);
    double alpha = plane_alpha (ff[], n);
      
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
static void embed_fraction_refine_css_test2 (Point point, scalar css_test2)
//void embed_fraction_refine_css_test2 (Point point, scalar css_test2) //2
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

    // coord n = facet_normal (point, css_test2, fss_test2);
    coord n = mycs (point, ff_oppo);
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
static void embed_face_fraction_fss_test2_refine_x (Point point, scalar s)
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

        // coord n = facet_normal (point, css_test2, fss_test2);
    coord n = mycs (point, ff_oppo);
    // double alpha = plane_alpha (css_test2[], n);
    double alpha = plane_alpha (ff_oppo[], n);
      
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





void restriction_Tlff (Point point, scalar s)
{
  double sum = 0.;
  double weight= 0.;
  foreach_child(){
    // if(ff[]>0.5){
    if(ff[]>0.0){
        // double local_weight = ff[];
        // sum += local_weight*s[];
        // weight += local_weight;
        sum += s[] ;
    }    
  }
  // s[] = sum/(1 << dimension)/(cm[] + 1e-30);
  if(weight>1e-12){
    // s[] = sum/weight;
    s[] = sum/(1 << dimension);
  }
}

void restriction_Tgff (Point point, scalar s)
{
  double sum = 0.;
  double weight= 0.;
  foreach_child(){
    // if(1.0-ff[]>0.5){
      if(1.0-ff[]>0.0){
        // double local_weight = 1.0-ff[];
        // sum += local_weight*s[];
        // weight += local_weight;
        sum += s[] ;
    }    
  }
  // s[] = sum/(1 << dimension)/(cm[] + 1e-30);
  if(weight>1e-12){
    // s[] = sum/weight;
    s[] = sum/(1 << dimension);
  }
}


///////////////////////////////////////////
/////////////////////////////////////
//refine function for css_test based scalar

void refine_embed_linear_css_test (Point point, scalar s)
//void refine_embed_linear2 (Point point, scalar s) //1
{
  foreach_child() {
    if (!css_test[]){
       s[] = Tsat00;
    }else {
      assert (coarse(css_test));
      int i = (child.x + 1)/2, j = (child.y + 1)/2;
#if dimension == 2
      if (coarse(fss_test.x,i) && coarse(fss_test.y,0,j) &&
	  (coarse(css_test) == 1. || coarse(css_test,child.x) == 1. ||
	   coarse(css_test,0,child.y) == 1. || coarse(css_test,child.x,child.y) == 1.)) {
	assert (coarse(css_test,child.x) && coarse(css_test,0,child.y));
	if (coarse(fss_test.x,i,child.y) && coarse(fss_test.y,child.x,j)) {
	  // bilinear interpolation
	  assert (coarse(css_test,child.x,child.y));
	  s[] = (9.*coarse(s) + 
		 3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
		 coarse(s,child.x,child.y))/16.;
	}
	else
	  // triangular interpolation	  
	  s[] = (2.*coarse(s) + coarse(s,child.x) + coarse(s,0,child.y))/4.;
      }
      else if (coarse(css_test,child.x,child.y) &&
	       ((coarse(fss_test.x,i) && coarse(fss_test.y,child.x,j)) ||
		(coarse(fss_test.y,0,j) && coarse(fss_test.x,i,child.y)))) {
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y))/4.;
      }
#else // dimension == 3
      int k = (child.z + 1)/2;
      if (coarse(fss_test.x,i) > 0.25 && coarse(fss_test.y,0,j) > 0.25 &&
	  coarse(fss_test.z,0,0,k) > 0.25 &&
	  (coarse(css_test) == 1. || coarse(css_test,child.x) == 1. ||
	   coarse(css_test,0,child.y) == 1. || coarse(css_test,child.x,child.y) == 1. ||
	   coarse(css_test,0,0,child.z) == 1. || coarse(css_test,child.x,0,child.z) == 1. ||
	   coarse(css_test,0,child.y,child.z) == 1. ||
	   coarse(css_test,child.x,child.y,child.z) == 1.)) {
	assert (coarse(css_test,child.x) && coarse(css_test,0,child.y) &&
		coarse(css_test,0,0,child.z));
	if (coarse(fss_test.x,i,child.y) && coarse(fss_test.y,child.x,j) &&
	    coarse(fss_test.z,child.x,child.y,k) &&
	    coarse(fss_test.z,child.x,0,k) && coarse(fss_test.z,0,child.y,k)) {
	  assert (coarse(css_test,child.x,child.y) && coarse(css_test,child.x,0,child.z) &&
		  coarse(css_test,0,child.y,child.z) &&
		  coarse(css_test,child.x,child.y,child.z));
	  // bilinear interpolation
	  s[] = (27.*coarse(s) + 
		 9.*(coarse(s,child.x) + coarse(s,0,child.y) +
		     coarse(s,0,0,child.z)) + 
		 3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
		     coarse(s,0,child.y,child.z)) + 
		 coarse(s,child.x,child.y,child.z))/64.;
	}
	else
	  // tetrahedral interpolation
	  s[] = (coarse(s) + coarse(s,child.x) + coarse(s,0,child.y) +
		 coarse(s,0,0,child.z))/4.;
      }
      else if (coarse(css_test,child.x,child.y,child.z) &&
	       ((coarse(fss_test.z,child.x,child.y,k) &&
		 ((coarse(fss_test.x,i) && coarse(fss_test.y,child.x,j)) ||
		  (coarse(fss_test.y,0,j) && coarse(fss_test.x,i,child.y))))
		||
		(coarse(fss_test.z,0,0,k) &&
		 ((coarse(fss_test.x,i,0,child.z) && coarse(fss_test.y,child.x,j,child.z)) ||
		  (coarse(fss_test.y,0,j,child.z) && coarse(fss_test.x,i,child.y,child.z))))
		||
		(coarse(fss_test.z,child.x,0,k) &&
		 coarse(fss_test.x,i) && coarse(fss_test.y,child.x,j,child.z))
		||
		(coarse(fss_test.z,0,child.y,k) &&
		 coarse(fss_test.y,0,j) && coarse(fss_test.x,i,child.y,child.z))
		))
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y,child.z))/4.;
#endif // dimension == 3
      else {
	// Pathological cases, use 1D gradients.
	s[] = coarse(s);
	foreach_dimension() {
	  if (coarse(fss_test.x,(child.x + 1)/2) && coarse(css_test,child.x))
	    s[] += (coarse(s,child.x) - coarse(s))/4.;
	  else if (coarse(fss_test.x,(- child.x + 1)/2) && coarse(css_test,- child.x))
	    s[] -= (coarse(s,- child.x) - coarse(s))/4.;
	}
      }
    }
  } 
}



/////////////////////////////////////////
////////////////////refine for scalar based on css_test2
//////////////////////////////////////

void refine_embed_linear_css_test2 (Point point, scalar s)
//void refine_embed_linear2 (Point point, scalar s) //1
{
  foreach_child() {
    if (!css_test2[]){
       s[] = Tsat00;
    }else {
      assert (coarse(css_test2));
      int i = (child.x + 1)/2, j = (child.y + 1)/2;
#if dimension == 2
      if (coarse(fss_test2.x,i) && coarse(fss_test2.y,0,j) &&
	  (coarse(css_test2) == 1. || coarse(css_test2,child.x) == 1. ||
	   coarse(css_test2,0,child.y) == 1. || coarse(css_test2,child.x,child.y) == 1.)) {
	assert (coarse(css_test2,child.x) && coarse(css_test2,0,child.y));
	if (coarse(fss_test2.x,i,child.y) && coarse(fss_test2.y,child.x,j)) {
	  // bilinear interpolation
	  assert (coarse(css_test2,child.x,child.y));
	  s[] = (9.*coarse(s) + 
		 3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
		 coarse(s,child.x,child.y))/16.;
	}
	else
	  // triangular interpolation	  
	  s[] = (2.*coarse(s) + coarse(s,child.x) + coarse(s,0,child.y))/4.;
      }
      else if (coarse(css_test2,child.x,child.y) &&
	       ((coarse(fss_test2.x,i) && coarse(fss_test2.y,child.x,j)) ||
		(coarse(fss_test2.y,0,j) && coarse(fss_test2.x,i,child.y)))) {
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y))/4.;
      }
#else // dimension == 3
      int k = (child.z + 1)/2;
      if (coarse(fss_test2.x,i) > 0.25 && coarse(fss_test2.y,0,j) > 0.25 &&
	  coarse(fss_test2.z,0,0,k) > 0.25 &&
	  (coarse(css_test2) == 1. || coarse(css_test2,child.x) == 1. ||
	   coarse(css_test2,0,child.y) == 1. || coarse(css_test2,child.x,child.y) == 1. ||
	   coarse(css_test2,0,0,child.z) == 1. || coarse(css_test2,child.x,0,child.z) == 1. ||
	   coarse(css_test2,0,child.y,child.z) == 1. ||
	   coarse(css_test2,child.x,child.y,child.z) == 1.)) {
	assert (coarse(css_test2,child.x) && coarse(css_test2,0,child.y) &&
		coarse(css_test2,0,0,child.z));
	if (coarse(fss_test2.x,i,child.y) && coarse(fss_test2.y,child.x,j) &&
	    coarse(fss_test2.z,child.x,child.y,k) &&
	    coarse(fss_test2.z,child.x,0,k) && coarse(fss_test2.z,0,child.y,k)) {
	  assert (coarse(css_test2,child.x,child.y) && coarse(css_test2,child.x,0,child.z) &&
		  coarse(css_test2,0,child.y,child.z) &&
		  coarse(css_test2,child.x,child.y,child.z));
	  // bilinear interpolation
	  s[] = (27.*coarse(s) + 
		 9.*(coarse(s,child.x) + coarse(s,0,child.y) +
		     coarse(s,0,0,child.z)) + 
		 3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
		     coarse(s,0,child.y,child.z)) + 
		 coarse(s,child.x,child.y,child.z))/64.;
	}
	else
	  // tetrahedral interpolation
	  s[] = (coarse(s) + coarse(s,child.x) + coarse(s,0,child.y) +
		 coarse(s,0,0,child.z))/4.;
      }
      else if (coarse(css_test2,child.x,child.y,child.z) &&
	       ((coarse(fss_test2.z,child.x,child.y,k) &&
		 ((coarse(fss_test2.x,i) && coarse(fss_test2.y,child.x,j)) ||
		  (coarse(fss_test2.y,0,j) && coarse(fss_test2.x,i,child.y))))
		||
		(coarse(fss_test2.z,0,0,k) &&
		 ((coarse(fss_test2.x,i,0,child.z) && coarse(fss_test2.y,child.x,j,child.z)) ||
		  (coarse(fss_test2.y,0,j,child.z) && coarse(fss_test2.x,i,child.y,child.z))))
		||
		(coarse(fss_test2.z,child.x,0,k) &&
		 coarse(fss_test2.x,i) && coarse(fss_test2.y,child.x,j,child.z))
		||
		(coarse(fss_test2.z,0,child.y,k) &&
		 coarse(fss_test2.y,0,j) && coarse(fss_test2.x,i,child.y,child.z))
		))
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y,child.z))/4.;
#endif // dimension == 3
      else {
	// Pathological cases, use 1D gradients.
	s[] = coarse(s);
	foreach_dimension() {
	  if (coarse(fss_test2.x,(child.x + 1)/2) && coarse(css_test2,child.x))
	    s[] += (coarse(s,child.x) - coarse(s))/4.;
	  else if (coarse(fss_test2.x,(- child.x + 1)/2) && coarse(css_test2,- child.x))
	    s[] -= (coarse(s,- child.x) - coarse(s))/4.;
	}
      }
    }
  } 
}


///////////////////////////////////
//////////restriction for scalar based on css_test
/////////////////////////////////

void restriction_embed_linear_css_test (Point point, scalar s)
//void restriction_embed_linear2 (Point point, scalar s)
{  
  // 0 children
  if (!css_test[] ) {
       s[] = Tsat00;
    return;
  }

  /**
  We first try to interpolate "diagonally". If enough child cells are
  defined (i.e. have non-zero embedded fractions), we return the
  corresponding value. */

  double val = 0., nv = 0.;
  for (int i = 0; i <= 1; i++)
#if dimension > 2
    for (int j = 0; j <= 1; j++)
#endif
      if (fine(css_test,0,i,j) && fine(css_test,1,!i,!j))
	val += (fine(s,0,i,j) + fine(s,1,!i,!j))/2., nv++;
  if (nv > 0.) {
    s[] = val/nv;
    return;
  }

  /**
  Otherwise, we use the average of the child cells which are defined
  (there is at least one). */
  
  coord p = {0.,0.,0.};
  foreach_child()
    if (css_test[])
      p.x += x, p.y += y, p.z += z, val += s[], nv++;
  assert (nv > 0.);
  s[] = val/nv;
}


///////////////////////////////////
//////////restriction for scalar based on css_test2
/////////////////////////////////

void restriction_embed_linear_css_test2 (Point point, scalar s)
//void restriction_embed_linear2 (Point point, scalar s)
{  
  // 0 children
  if (!css_test2[] ) {
       s[] = Tsat00;
    return;
  }

  /**
  We first try to interpolate "diagonally". If enough child cells are
  defined (i.e. have non-zero embedded fractions), we return the
  corresponding value. */

  double val = 0., nv = 0.;
  for (int i = 0; i <= 1; i++)
#if dimension > 2
    for (int j = 0; j <= 1; j++)
#endif
      if (fine(css_test2,0,i,j) && fine(css_test2,1,!i,!j))
	val += (fine(s,0,i,j) + fine(s,1,!i,!j))/2., nv++;
  if (nv > 0.) {
    s[] = val/nv;
    return;
  }

  /**
  Otherwise, we use the average of the child cells which are defined
  (there is at least one). */
  
  coord p = {0.,0.,0.};
  foreach_child()
    if (css_test2[])
      p.x += x, p.y += y, p.z += z, val += s[], nv++;
  assert (nv > 0.);
  s[] = val/nv;
}

void conservative_refine (Point point, scalar s)
{  
  double cc = ff[];
  double scc = s[];
  if (cc <= 0. || cc >= 1.)
    refine_bilinear(point,s);
  else {

    coord n = mycs (point, ff);
    double alpha = plane_alpha (cc, n);

    coord a, b;
    foreach_dimension() {
      a.x = 0.; b.x = 0.5;
    }   
        
    foreach_child() {
      coord nc; 
      foreach_dimension()
        nc.x = child.x*n.x;
      double crefine = rectangle_fraction (nc, alpha, a, b); 
      if (s.inverse)
        s[] = scc/(1. - cc)*(1. - crefine);
      else
        s[] = scc/cc*crefine;
    }   
  }   
}

//aiml aimg refine
void refine_aim(Point point, scalar s){
  double cc = css_test3[];
  double s_value = s[];
    if(cc<0 && cc>1.0){
          foreach_child(){
              s[] = 0;
          }
    }else{
        coord n = mycs (point, css_test3);
        double alpha = plane_alpha (cc, n);
        //for alpha_refine
        coord m;
        double alphac = 2.*alpha;
        foreach_dimension()
              m.x = n.x;

        foreach_child(){
              coord p_c_child_c; // tansfer to parent's coordate
              double cc2;
              double alpha_child;

              static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
              coord nc;
              foreach_dimension(){
                  nc.x = child.x * n.x;
              }
              cc2 = rectangle_fraction(nc,alpha,a,b);  //child_ff3
              alpha_child=alphac;
              if(cc2>0.0 && cc2<1.0){
                  s[] = s_value;  
              }else{
                  s[] = 0;
              }
        }
    }
}

void restriction_aim(Point point, scalar s){
    double volume_total=0.0;
    foreach_child(){
        volume_total += css_test3[];
    }
    volume_total = volume_total /(1<<dimension);
    if(volume_total>0.0 && volume_total<1.0){
        double weight=0.0;
        double value_tot=0.0;
        foreach_child(){
            if(css_test3[]>0.0 && css_test3[]<1.0){
                double weight_local = 1;
                weight += weight_local;
                value_tot += (weight_local*s[]);
            }
        }
        if(weight>0){
            s[] = value_tot/weight;
        }
    }
}

void across_interface(scalar f_old, scalar f_new, scalar TT){
  foreach(){
    if(css_test3[]<=0.0){
        if((f_old[]<0.5 && f_new[]>=0.5) || (f_old[]>=0.5 && f_new[]<0.5)){
            TT[] = Tsat00;
        }
    }
  }

}



struct smooth1 {
  scalar df;
  scalar tempf;
};
void smooth_for_arm(struct smooth1 q){
  scalar smf[];
   scalar ssmf[];
   scalar sssmf[];
   scalar ssssmf = q.df;
   scalar tempf = q.tempf;
  // scalar ssssmf[];
  // scalar sssssmf[];
   //for (int ii=1;ii<=2;ii++) {

#if dimension == 2
    foreach()
        smf[] = (4.*tempf[] +
          2.*(tempf[0,1] + tempf[0,-1] + tempf[1,0] + tempf[-1,0]) +
          tempf[-1,-1] + tempf[1,-1] + tempf[1,1] + tempf[-1,1])/16.;


    foreach()
        ssmf[] = (4.*smf[] +
          2.*(smf[0,1] + smf[0,-1] + smf[1,0] + smf[-1,0]) +
          smf[-1,-1] + smf[1,-1] + smf[1,1] + smf[-1,1])/16.;

    foreach()
        sssmf[] = (4.*ssmf[] +
          2.*(ssmf[0,1] + ssmf[0,-1] + ssmf[1,0] + ssmf[-1,0]) +
          ssmf[-1,-1] + ssmf[1,-1] + ssmf[1,1] + ssmf[-1,1])/16.;

    foreach()
        ssssmf[] = (4.*sssmf[] +
          2.*(sssmf[0,1] + sssmf[0,-1] + sssmf[1,0] + sssmf[-1,0]) +
          sssmf[-1,-1] + sssmf[1,-1] + sssmf[1,1] + sssmf[-1,1])/16.;

    // foreach()
    //     sssssmf[] = (4.*ssssmf[] +
    //       2.*(ssssmf[0,1] + ssssmf[0,-1] + ssssmf[1,0] + ssssmf[-1,0]) +
    //       ssssmf[-1,-1] + ssssmf[1,-1] + ssssmf[1,1] + ssssmf[-1,1])/16.;

#else  //dimension == 3
     foreach(){
          double ccc[3][3][3];
          for(int ii=-1;ii<=1;ii++){
              for(int jj=-1;jj<=1;jj++){
                  for(int kk=-1;kk<=1;kk++){
                       ccc[ii+1][jj+1][kk+1]=tempf[ii,jj,kk];
                  }//kk
              } //jj
          } //ii
          smf[]=smff(ccc);
     }
   //  smf.refine = refine_injection;
   //  smf.restriction=restriction_average;
   //  boundary({smf});
     printf("smff 1 pid=%d\n",pid());
     foreach(){
         // printf("pid=%d:%g %g %g     ",pid(),x,y,z);
          double ccc[3][3][3];
          for(int ii=-1;ii<=1;ii++){
              for(int jj=-1;jj<=1;jj++){
                  for(int kk=-1;kk<=1;kk++){
                       ccc[ii+1][jj+1][kk+1]=smf[ii,jj,kk];
                  }//kk
              } //jj
          } //ii
          ssmf[]=smff(ccc);
          //printf("finish is_local=%d,pid=%d:%g %g %g\n",is_local(cell),pid(),x,y,z);
     }
      printf("smff 2 pid=%d\n",pid());
     foreach(){
          double ccc[3][3][3];
          for(int ii=-1;ii<=1;ii++){
              for(int jj=-1;jj<=1;jj++){
                  for(int kk=-1;kk<=1;kk++){
                       ccc[ii+1][jj+1][kk+1]=ssmf[ii,jj,kk];
                  }//kk
              } //jj
          } //ii
          sssmf[]=smff(ccc);
     }
      printf("smff 3 pid=%d\n",pid());
     foreach(){
          double ccc[3][3][3];
          for(int ii=-1;ii<=1;ii++){
              for(int jj=-1;jj<=1;jj++){
                  for(int kk=-1;kk<=1;kk++){
                       ccc[ii+1][jj+1][kk+1]=sssmf[ii,jj,kk];
                  }//kk
              } //jj
          } //ii
          ssssmf[]=smff(ccc);
     }
      printf("smff 4 pid=%d\n",pid());

#endif

}

// #include <stdio.h>
// #include <math.h>

// Normalize a vector
void normalizeVector(double* vx, double* vy) {
    double magnitude = sqrt((*vx) * (*vx) + (*vy) * (*vy));
    if (magnitude != 0) {
        *vx /= magnitude;
        *vy /= magnitude;
    }
}

// Function to transform a vector into a new coordinate system defined by the normal
void transformVectorToSurfaceCoordinate(double ux, double uy, double nx, double ny, double* un, double* ut) {
    // Ensure the normal vector is normalized
    normalizeVector(&nx, &ny);

    // Calculate un and ut without modifying ux, uy
    *un = ux * nx + uy * ny;
    *ut = ux * (-ny) + uy * nx;
}

// Function to transform from surface normal coordinates back to x-y coordinates
void transformSurfaceCoordinateToVector(double un, double ut, double nx, double ny, double* ux, double* uy) {
    // Ensure the normal vector is normalized
    normalizeVector(&nx, &ny);

    // Calculate ux and uy from un and ut
    *ux = un * nx + ut * (-ny);
    *uy = un * ny + ut * nx;
}




// extern scalar topo_mask_s;
extern int maxl;
// scalar u, scalar css_test3_n, solid_m
void get_u_ghost(vector u){
      vector u_surface[];
      // double lambdas = (L0)/(1<<maxl);
      double lambdas = HUGE; //slip boundary

      scalar css_test3_n_neg[];
      foreach(){
        css_test3_n_neg[] = 1.0-css_test3_n[];
      }
      scalar solid_neg_alpha[];
      vector solid_n[];
      reconstruction (css_test3_n_neg, solid_n, solid_neg_alpha);


      // u -> u_surface
      foreach(){
         if(topo_mask_s[]==0){
              double ux = u.x[]; // Vector u in x-y coordinates
              double uy = u.y[];
              double nx = solid_n.x[]; // Updated normal of the solid surface
              double ny = solid_n.y[];
              double un, ut; // Components of u in the surface normal coordinate system

              // Transform u to the surface normal coordinate system
              transformVectorToSurfaceCoordinate(ux, uy, nx, ny, &un, &ut);
              u_surface.x[] = un; //
              u_surface.y[] = ut;
         }
      }

  //     *n = facet_normal (point, css_test3_n, fs);
  // double alpha = plane_alpha (css_test3_n[], *n);
  // double area = plane_area_center (*n, alpha, p);
  // normalize (n);
  // return area;

      vector o3[],o1[];
      foreach(){ //get o3, center of interface
        coord n_css_test3_n,o3_temp;
        foreach_dimension(){
          o3.x[]=0.0;
          n_css_test3_n.x = -solid_n.x[];
          o3_temp.x = 0.0;
        }
        if(css_test3_n[]>0.0 && css_test3_n[]<1.0){
            double alpha = plane_alpha (css_test3_n[], n_css_test3_n);
            double area = plane_area_center (n_css_test3_n, alpha, &o3_temp);
            o3.x[]=x + o3_temp.x*Delta;
            o3.y[]=y + o3_temp.y*Delta;
        }
      }

    foreach(){ //get o1center of the fluid volume
      coord n_css_test3_n,o1_temp;
      foreach_dimension(){
        o1.x[]=0.0;
        n_css_test3_n.x = -solid_n.x[];
        o1_temp.x = 0.0;
      }
      if(css_test3_n[]>0.0 && css_test3_n[]<1.0){
         double alpha = plane_alpha (css_test3_n[], n_css_test3_n);
         plane_center (n_css_test3_n, alpha, css_test3_n[], &o1_temp);
         o1.x[]=x + o1_temp.x*Delta;
         o1.y[]=y + o1_temp.y*Delta;
      }
    }

    // get a1
    scalar a1[];
    foreach(){
        a1[] = 0;
        // project
        if(css_test3_n[]>0.0 && css_test3_n[]<1.0){
            coord o1_o3;
            foreach_dimension(){
               o1_o3.x = o1.x[] - o3.x[]; 
            }
            double n_mod;
            n_mod = sqrt(sq(solid_n.x[])+sq(solid_n.y[]));
            if(n_mod>1e-20){
              double dott = (o1_o3.x*solid_n.x[] + o1_o3.y*solid_n.y[])/n_mod;
              a1[] = fabs(dott);
            }
        }
    }

    scalar solid_u[];
    foreach(){
      if(topo_mask_s[]==0){
        solid_u[] = u_surface.y[]*a1[]/(a1[]+lambdas); //t direction
      }
    }

    // for topo_mask_s==1
    foreach(){
      if(topo_mask_s[]==1 && level==level_interface){ //
          double totalx =0.0;
          double weightx = 0.0;

          double totaly =0.0;
          double weighty = 0.0;
         //x minus side
         //x plus side
         if(topo_mask_s[1,0]==0){
              coord u_surface_local;
              //get a2
              //solid_n.x*x+solid_n.y*y=alpha, origin is x,y of   Point[1,0]
              // current point coordinate then is (-1,0)
              double n_mod = sqrt(sq(solid_n.x[1,0])+sq(solid_n.y[1,0]));
              if(n_mod>1e-20){
                double a2 = fabs(solid_n.x[1,0]*(-1)+solid_n.y[1,0]*0-solid_neg_alpha[1,0])/n_mod*Delta;
                u_surface_local.x = -(a2)/(a1[1,0])*u_surface.x[1,0];//n
                u_surface_local.y = u_surface.y[1,0]*(lambdas-a2)/(a1[1,0]+lambdas);//
                double ux; // Vector u in x-y coordinates
                double uy;
                double nx = solid_n.x[1,0]; // Updated normal of the solid surface
                double ny = solid_n.y[1,0];
                double un = u_surface_local.x;
                double ut = u_surface_local.y; // Components of u in the surface normal coordinate system

                transformSurfaceCoordinateToVector(un, ut, nx, ny, &ux, &uy);
              
                totalx += ux * 1;
                weightx += 1;

                totaly += uy * 1;
                weighty += 1;
              }
         }
         //y minus side
         //y plus side
        
        if(weightx>0.0 && weighty>0.0){
                u.x[] = totalx/weightx;
                u.y[] = totaly/weighty;
         }
      }else if(topo_mask_s[]>=2){
          u.x[] = 0;
          u.y[] = 0;
      }
    }


}




void get_ps_ghost(scalar ps){
      scalar p_surface[];

      scalar css_test3_n_neg[];
      foreach(){
        css_test3_n_neg[] = 1.0-css_test3_n[];
      }
      scalar solid_neg_alpha[];
      vector solid_n[];
      reconstruction (css_test3_n_neg, solid_n, solid_neg_alpha);


      // // u -> u_surface
      // foreach(){
      //    if(topo_mask_s[]==0){
      //         double ux = u.x[]; // Vector u in x-y coordinates
      //         double uy = u.y[];
      //         double nx = solid_n.x[]; // Updated normal of the solid surface
      //         double ny = solid_n.y[];
      //         double un, ut; // Components of u in the surface normal coordinate system

      //         // Transform u to the surface normal coordinate system
      //         transformVectorToSurfaceCoordinate(ux, uy, nx, ny, &un, &ut);
      //         u_surface.x[] = un; //
      //         u_surface.y[] = ut;
      //    }
      // }

  //     *n = facet_normal (point, css_test3_n, fs);
  // double alpha = plane_alpha (css_test3_n[], *n);
  // double area = plane_area_center (*n, alpha, p);
  // normalize (n);
  // return area;

      vector o3[],o1[];
      foreach(){ //get o3, center of interface
        coord n_css_test3_n,o3_temp;
        foreach_dimension(){
          o3.x[]=0.0;
          n_css_test3_n.x = -solid_n.x[];
          o3_temp.x = 0.0;
        }
        if(css_test3_n[]>0.0 && css_test3_n[]<1.0){
            double alpha = plane_alpha (css_test3_n[], n_css_test3_n);
            double area = plane_area_center (n_css_test3_n, alpha, &o3_temp);
            o3.x[]=x + o3_temp.x*Delta;
            o3.y[]=y + o3_temp.y*Delta;
        }
      }

    foreach(){ //get o1center of the fluid volume
      coord n_css_test3_n,o1_temp;
      foreach_dimension(){
        o1.x[]=0.0;
        n_css_test3_n.x = -solid_n.x[];
        o1_temp.x = 0.0;
      }
      if(css_test3_n[]>0.0 && css_test3_n[]<1.0){
         double alpha = plane_alpha (css_test3_n[], n_css_test3_n);
         plane_center (n_css_test3_n, alpha, css_test3_n[], &o1_temp);
         o1.x[]=x + o1_temp.x*Delta;
         o1.y[]=y + o1_temp.y*Delta;
      }
    }

    // get a1
    scalar a1[];
    foreach(){
        a1[] = 0;
        // project
        if(css_test3_n[]>0.0 && css_test3_n[]<1.0){
            coord o1_o3;
            foreach_dimension(){
               o1_o3.x = o1.x[] - o3.x[]; 
            }
            double n_mod;
            n_mod = sqrt(sq(solid_n.x[])+sq(solid_n.y[]));
            if(n_mod>1e-20){
              double dott = (o1_o3.x*solid_n.x[] + o1_o3.y*solid_n.y[])/n_mod;
              a1[] = fabs(dott);
            }
        }
    }


    // for topo_mask_s==1
    foreach(){
      if(topo_mask_s[]==1 && level==level_interface && topo_mask[]<=0){ //
          double totalx =0.0;
          double weightx = 0.0;

          bool flag=false;
          foreach_neighbor(){
            if(intersect_true[]==1){
                flag=true;
            }
          }
         //x minus side
         //x plus side
         if(topo_mask_s[1,0]==0 && flag){ // (p[1]-p[])/delta/rho*dt=u(t+1)-u(t)
              double p_surface_local;
              //get a2
              //solid_n.x*x+solid_n.y*y=alpha, origin is x,y of   Point[1,0]
              // current point coordinate then is (-1,0)
              double n_mod = sqrt(sq(solid_n.x[1,0])+sq(solid_n.y[1,0]));
              if(n_mod>1e-20){
                double a2 = fabs(solid_n.x[1,0]*(-1)+solid_n.y[1,0]*0-solid_neg_alpha[1,0])/n_mod*Delta;
                p_surface_local = -(a2)/(a1[1,0])*ps[1,0];//n
                totalx += p_surface_local * 1;
                weightx += 1;
              }
         }
         //y minus side
         //y plus side
        
        if(weightx>0.0){
                ps[] = totalx/weightx;
         }
      }
    }


}
