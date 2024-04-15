extern scalar ff;
extern int globali;
extern int minl;
//extern scalar cs;
//extern face vector fs;
extern scalar res_ps;
extern scalar topo_mask_s;

extern bool flag_topos_advect_uf;
extern scalar intersect_true;
struct Linear_ps {
  face vector alpha_temp;
  scalar ps;
  face vector usf;
  scalar topo_mask;
  scalar source;
  int level_interface;   
  double dt; 
};

void poisson_ps_usf1 (struct Linear_ps q) {
   
        //rho sv
     //set scalar A3,A4, vector Alef, Arigh,
     face vector alpha_temp = q.alpha_temp;
     scalar ps = q.ps;
     face vector usf = q.usf;
     scalar topo_mask = q.topo_mask;
     scalar source = q.source;
     int level_interface = q.level_interface; 
     double dt = q.dt;
 
/*
    if(1==1){
      char name82[80];  
      sprintf(name82,"pid%d-i%d",pid(),globali);
      FILE * fp11= fopen(name82,"w");
      foreach_level(7){
          fprintf(fp11,"%g %g %g\n",x,y,z);
      }
      fclose(fp11);
    }
*/
     //printf("pid=%d begin usf1-lin34\n",pid());


     //it seems that point could not be used in foreach and foreach_dimension

//printf("pid=%d begin usf1-line74 depth()=%d\n",pid(),depth());
     vector Alef[],Arigh[];
     scalar A3[],A4[];
     foreach(){
        foreach_dimension(){
           Alef.x[]=0.0;
           Arigh.x[]=0.0;
        }
        A3[]=0.0;
        A4[]=0.0;
     }
     //for(int l=depth();l>=1;l--){
     //    boundary_level({A3,A4},Alef.x,Arigh.x,l);
     //}
 //   printf("pid=%d begin usf1-line89 depth()=%d\n",pid(),depth()); 
/*
     for(int l=depth();l>=1;l--){
          foreach_halo(prolongation,l-1){
              Alef.x[] = 0.0;
              Arigh.x[]= 0.0;
              A3[]=0.0;
              A4[]=0.0;
          
          Point neip = neighborp(1);
          if(is_boundary_box2(neighbor(1))){
              foreach_neighbor(1){
                   if(point.i==neip.i && point.j==neip.j && point.k==neip.k){
                        A3[] = 0.0;
                        A4[] = 0.0;
                        Alef.x[]=0.0;
                        Arigh.x[]=0.0;
                   }
              }

            }
          }
          boundary_level({A3,A4,Alef.x,Arigh.x},l);
     }
 */
 //printf("pid=%d begin usf1-line112 depth()=%d\n",pid(),depth());

     int numbercells=0;
  //   bool get_delta_interface=false;
     double Delta_interface = HUGE;
     scalar numbercells3[];
     foreach(){
        numbercells3[]=0;
     }
     foreach(reduction(+:numbercells)){
      //   if(level == level_interface && (cs[]>0)){
          if(level == level_interface && (cs[]>0)){
	     double temp_A3=0.0; 
	     int current = topo_mask[];
	         if(current==-1 || current==-2){
                       numbercells+=1;
                       numbercells3[]=1;
                       bool lefedge= true; //(fs.x[]>0 && fs.x[]<=1); // if solid(false), then wall boundary
                       bool righedge= true; //(fs.x[1]>0 && fs.x[1]<=1);
                       foreach_dimension(){
                              //Alef.x[]= 2.0*facemask.x[]/(Delta*Delta)/(rhov[]+rhov[-1]);
                              //Arigh.x[]= 2.0*facemask.x[1]/(Delta*Delta)/(rhov[1]+rhov[]);
                          if(lefedge){
                              Alef.x[]= dt/(Delta*Delta)*alpha_temp.x[];
                          }
                          if(righedge){
                              Arigh.x[]= dt/(Delta*Delta)*alpha_temp.x[1];
                          }
                             temp_A3 += Alef.x[]+Arigh.x[];
                       }
                       A3[] = temp_A3;
	               //A4[]=  -source_pc[]*cm[]; //times cm[] for axi-symetric z 
                  //  A4[]=  -source[];//*cm[]; //20230619
                 A4[]= 0.0;     //20230928
             }else if(current==0){
                      numbercells+=1;
                      numbercells3[]=1;
                      foreach_dimension(){
                            int lef=topo_mask[-1];
                            int righ=topo_mask[1];
                            bool lefedge= true; //(fs.x[]>0 && fs.x[]<=1); // if solid(false), then wall boundary
                            bool righedge= true;// (fs.x[1]>0 && fs.x[1]<=1);
                              if((lef==0 || lef==-1) && lefedge){
                              //Alef.x[]= 2.0*facemask.x[]/(Delta*Delta)/(rhov[]+rhov[-1]);
                                                   Alef.x[]= dt/(Delta*Delta)*alpha_temp.x[];
                              }
                              if((righ==0 || righ==-1) && righedge){
                              //Arigh.x[]= 2.0*facemask.x[1]/(Delta*Delta)/(rhov[1]+rhov[]);
                                                   Arigh.x[]= dt/(Delta*Delta)*alpha_temp.x[1];
                              }
                            temp_A3 += Alef.x[]+Arigh.x[];
                      }  
                      A3[] = temp_A3;
	             // A4[]=  -source_pc[]*cm[]; //times cm[] for axi-symetric  
                  if(!flag_topos_advect_uf){
                     A4[]=  -source[];//*cm[];  //20230619
                  }else{
                     if(topo_mask_s[]==0 && (intersect_true[]==0)){
                     // if(topo_mask_s[]==0){
                        A4[]=  0;//*cm[];  //20230619
                     }else{
                        A4[]=  -source[];//*cm[];  //20230619
                     }
                  }
             }else{
               //p=0.0
                foreach_dimension(){
                     Alef.x[]=0.0;
                     Arigh.x[]=0.0;
                }
                A3[]=1.0;
                A4[]=0.0; 
            }
		      // if(!get_delta_interface){
            //      get_delta_interface = true;
            //      Delta_interface = Delta;
            // }
        } //level == level_interface

     } // foreach

//printf("pid=%d begin usf1-line173 depth()=%d\n",pid(),depth());

///////////////////// linear solver
     int maxitt=2000,beta1=1.2;
     //double ps_tolerance = 1.0/max(Trhol,Trhog)/(Delta_interface*Delta_interface)*1e-7*dt;
    // double ps_tolerance = 0.000000001/dt; //0.0001/dt;
    double ps_tolerance = 1e-4;//0.000001;//1e-4;//0.000001;//0.00000001;//0.00000001;//0.00000001/dt; //0.0001/dt;
     foreach(){
        ps[] = 0.0;
     }
     foreach_face(){
        usf.x[] = 0.0;
     }
     boundary({ps,usf});
//#define output_residual_ps 1;
#ifdef output_residual_ps
  char name65[80];
  sprintf(name65,"residual_ps_%g.dat", t);
  FILE * fp65 = fopen(name65,"a");
  fprintf (fp65, "# it-residual tres2 perror=tres2*maxrho*delta*delta\n");

#endif  
     int itt; 
    // clock_t start1,end1;
    // double during1;
    // start1=clock();
     for(itt=1;itt<=maxitt;itt++){
    //     if(itt%5==0){
    //        end1 = clock();
    //        during1 =  (end1 - start1)/(double)CLOCKS_PER_SEC;
    //        fprintf(stderr,"itt=%d elapse=%g\n",itt,during1);
    //     }
         foreach(){
            int temp3;
            temp3 = topo_mask[];
            //if( temp3==3 || temp3==-3){
            if( temp3>=1 || temp3==-3){
                ps[]=0.0;
            }
           
         }
        // boundary({ps});
         /*
         if(itt%5==0 && pid()==1){
            FILE * file11 = fopen("huihui.dat","a");
            fprintf(file11,"pid=%d itt=%d\n",pid(),itt);
            fclose(file11);
         }
         */
//20220814  add for periodic but not sure for other condition
         //boundary({ps}); 
         //boundary_level({ps},level_interface) 
         // for(int l=depth();l>=0;l--){
         //     boundary_level({ps},l);
         // }
         boundary({ps});
         // get_ps_ghost(ps);
         foreach(){   
            if((level==level_interface) && (numbercells3[]==1)){
               double part1=0.0;
               double part2;
               foreach_dimension(){
                     part1 = part1 +  (beta1/A3[])*(Alef.x[]*ps[-1] + Arigh.x[]*ps[1]);
               }  
               part2=(1.0-beta1)*ps[] + part1 + (beta1/A3[])*A4[];
               ps[] = part2;
            }
          }  
          boundary({ps});
         //  get_ps_ghost(ps);
         // boundary_level({ps},level_interface);

         double T_res2=0.0;
         int numbercells2=0;
         double maxres_ps=0.0;
         foreach(reduction(+:T_res2) reduction(+:numbercells2) reduction(max:maxres_ps)){
            res_ps[] = 0.0;
            if((level == level_interface) && (numbercells3[]==1)){
              numbercells2+=1;
              double temp_res2=0.0;
              //temp_res2 = -T[]*A3[] + A1[]*T[-1] + A2[]*T[1] + A4[];
              foreach_dimension(){
                  temp_res2 = temp_res2 + Alef.x[]*ps[-1] + Arigh.x[]*ps[1];
              }
              temp_res2 = temp_res2 -ps[]*A3[] + A4[];
              T_res2=T_res2+pow(temp_res2,2.0); 
              res_ps[] = temp_res2/dt;
              if(fabs(temp_res2)>maxres_ps){
                  maxres_ps = fabs(temp_res2);
              }
              
            }
         }
	   //  if(numbercells>=1){
      //        T_res2 = T_res2/numbercells;
      //   }else{
      //        T_res2 = 0.0;
      //   }
      if(numbercells2<=0){
            printf("numbercell2=%d numbercells=%d\n",numbercells2,numbercells);
      }
      T_res2 = T_res2/numbercells2;
          // catch divergence
         int diverged=0;
         foreach(reduction(max:diverged)){
            if(level == level_interface){
             foreach_dimension(){
                if( (Alef.x[] != Alef.x[]) || 
                  (Arigh.x[] != Arigh.x[])) {
                    diverged=1;
                }
             }
             if((diverged==1) || (ps[]!=ps[])){
  //          char name5[80];
  //          sprintf(name5,"poisson_coeff_%g.dat",t);
  //          FILE * fp5 = fopen(name5,"a");
  //          fprintf (fp5, "I=%d J=%d ",point.i,point.j);
  //          fprintf (fp5, "A or p is NaN after %d iterations \n",itt);
  //          foreach_dimension()
  //             fprintf (fp5, "A1=%g A2=%g ", lefA.x[], righA.x[]);
  //          fprintf (fp5, "A3=%g A4=%g T=%g \n", AA3[], AA4[], T[]);
  //          fprintf (fp5, "A or p is NaN after %d iterations \n",itt);
  //          fclose(fp5);
               //  fprintf(stderr,"ps_poisson_error\n");
           
             }
            } //if level_interface

         }
          
   

          if(diverged==1){
            printf("diverged for ps; pid=%d\n",pid());
              exit(1);
           }

      //if ((T_res2*numbercells) > (1e16) ) {
      if ((T_res2*numbercells2) > (1e16) ) {
            //char name64[80];
            //sprintf(name64,"poisson_ps_%g.dat", t);
           // FILE * fp64 = fopen(name64,"a");
            
            //fprintf (fp64, "Poisson solver diverged after %d iterations \n",itt);
            fprintf (stderr, "Poisson solver diverged after %d iterations \n",itt);
            //fclose(fp64);

            // char name66[80];
            // sprintf(name66,"poisson_ps_A%g.dat", t);

            // FILE * fp66 = fopen(name66,"a");
            // foreach(){
            //   if(fabs(topo_mask[])<3)        
            //        fprintf(fp66,"%g %g %g %g %g %g %g %g %g %g\n",topo_mask[],x,y,z,Alef.x[],Arigh.x[],Alef.y[],Arigh.y[],A3[],A4[]);
            // }
            // fclose(fp66);


 //           exit(1);
       }
   
         double tres2=sqrt(T_res2);
	  // FILE * fp3 = fopen(name3,"a");
	 //double ittt=pow((double)(N),2.0);
#ifdef output_residual_ps
             double perror= tres2 *max(Trhol,Trhog)*(Delta*Delta);
             fprintf(fp65,"%d %g %g %d %d\n",itt,tres2,perror,diverged,numbercells);
#endif
	 //  fprintf(fp3,"%d %g %d ittt=%g \n",itt,tres2,diverged,ittt);
	 //  fclose(fp3);
	  if (tres2<ps_tolerance){
	       break;
	  }


         // // // if(pid()==0){

         // // //    char name66[80];
         // // //    sprintf(name66,"poisson_ps_A%g.dat", t);

         // // //    FILE * fp66 = fopen(name66,"a");
         // // //    foreach(){
         // // //      if(fabs(topo_mask[])<3)        
         // // //           fprintf(fp66,"%d %g %g\n",itt,maxres_ps/dt,tres2);
         // // //    }
         // // //    fclose(fp66);
         // // // }

         
         // if(maxres_ps/dt<ps_tolerance){
         //       break;
         // }

//20220814  add for periodic but not sure for other condition
      boundary({ps});
      // get_ps_ghost(ps);

      // if(pid()==0){
      //       printf("    itt=%d\n",itt);
      // }
     }//itt
fprintf(stderr,"temperature itt=%d\n",itt);
#ifdef output_residual_ps
     fclose(fp65);
#endif

    

  //update for usf
   foreach_face(){
      if(level == level_interface){
      //   if(facemask.x[]==1 && fs.x[]>0){
         if(fs.x[]>0){
            if(topo_mask[]==-1 || 
                (topo_mask[]==-2 && topo_mask[-1]==-1 ) || 
                (topo_mask[]==0  && ( topo_mask[-1]==0 || topo_mask[-1]==-1))){
                    usf.x[] = -dt*alpha_temp.x[]*face_gradient_x (ps, 0);
            }
        } 
      }// if level_interface
   }
   foreach_face(){
      if(fs.x[]==0.0){
          usf.x[] = 0.0;
      }
   }
}

#define EPS32 0.0

struct Topo_m2 {
    scalar topo_mask;
};
void get_topomask(struct Topo_m2 q) {
     scalar topo_mask = q.topo_mask;

     foreach(){
        int phase_sign = (ff[]>=0.5-EPS32) ? 1 : -1;
	topo_mask[] = 3*phase_sign;
        //if(ff[]<1.0-EPS32 && ff[]>EPS32){
         if(ff[]<1.0-EPS32 && ff[]>EPS32 && (css_test3_n[]>0.0)){ // topo_mask==0 could not be inside of the solid
         // if(ff[]<1.0-EPS32 && ff[]>EPS32){ //could be inside of the solid
 		topo_mask[] = 0;
        }
     }
     boundary({topo_mask});
     foreach(){
       // if(level==level_interface){
         // if(level==level_interface && (css_test3_n[]>0.0)){
         if(level==level_interface){
            bool is1= false;
            foreach_dimension(){
               int temp1=topo_mask[1];
               int temp2=topo_mask[-1];
            if(temp1==0 || temp2==0){
               is1 = true;
               }
            }
            int temp3=topo_mask[];
            if (is1 && temp3!=0){
              // topo_mask[] = (ff[]>=0.5-EPS32) ? 1 : -1 ;
	      if(temp3==3){
	          topo_mask[] = 1;
	      }else if(topo_mask[]==-3){
	          topo_mask[] = -1;
	      }
            }
        }
        
     } 
  for(int phase=0;phase<=1;phase++){
     foreach(){
        // if(level==level_interface){
         //  if(level==level_interface && (css_test3_n[]>0.0)){
         if(level==level_interface ){
               int phase_sign = ff[]>=0.5-EPS32 ? 1 : -1;
            bool is1= false;
               foreach_dimension(){
                  int temp1=topo_mask[1];
                  int temp2=topo_mask[-1];
               if( temp1==(2*phase-1) || temp2==(2*phase-1)){
                  is1 = true;
                  }
               }
               int temp3 = topo_mask[];
               if (is1 && temp3==3*phase_sign){
                  topo_mask[] = 2*(2*phase-1) ;
               }
         }
     }
  }


       
}
