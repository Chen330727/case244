//The flowing calculate temperature gradient at the interface using method in Paris
//this method could be used in (1) mass_transfer_rate when calculate teperature gradient at the interface
//which is related to source term (2) in residual3 of poisson3 when calculate diffustion flux 
//acrros the interface
//Note: ff5 should be the original ff; 
#include "fractions.h"
#include "curvature.h"  
#include "gradient_poisson.h"
//#include "geometry.h"
extern scalar topo_mask;
extern scalar topo_mask_s;
extern int level_interface,globali;
extern scalar css_test3_n,css_test3,css_test,css_test2;
extern scalar T;
extern double Tsat00;


extern scalar T_grad_method_g;
extern scalar T_grad_method_l;
extern scalar Tl,Tg;
extern scalar grad_embed_3rd_g,grad_embed_3rd_l;
//extern bool new_flag;

extern scalar intersect_true;


// Function prototypes
bool linear_regression(double *x, double *y, int n, double *m, double *b);
bool quadratic_regression(double *x, double *y, int n, double *a, double *b, double *c);
double predict_linear(double m, double b, double x);
double predict_quadratic(double a, double b, double c, double x);
bool are_points_too_close(double *distances, int n, double delta);

coord interface_normal7 (Point point, scalar c, vector h)
{
  coord n;
  if (!h.x.i || (n = height_normal (point, c, h)).x == nodata)
    n = mycs (point, c);
  return n;
}

//note::!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//phaseg, for gas T[]-Tsat, out of liquid(not out of gas)
//phasel, for liquid Tsat-T[], out of liquid 

struct Tgradleon {
  scalar ff5;
  scalar phaseg;
  scalar phasel;
  scalar aiml;
  scalar aimg;
};


// void Tgrad_leon(scalar ff5,scalar phaseg,scalar phasel){
// void Tgrad_leon(struct Tgradleon q){
void Tgrad_leon_comment(struct Tgradleon q){
  scalar ff5 = q.ff5;
  scalar phaseg = q.phaseg;
  scalar phasel = q.phasel;
//   double length_threshold = 0.05;
 double length_threshold = 0.01;
 bool is_constant_m=false;//;true;//false;//true;
    //ff5 should be original ff
 //       printf("pid=%d, T_leon begin\n",pid());
    vector hh2[];
 //   printf("pid=%d, T_leon 22\n",pid());
    heights(ff5,hh2);
 //   printf("pid=%d, T_leon 24\n",pid());
    scalar distanceT[];
    scalar get_novalue_flag[];
  //1: calculate Tgradient for pure cells
      //require topo_mask
    foreach(){
        phaseg[]=0.0;
        phasel[]=0.0;
        distanceT[] = HUGE;
    }
    double threshold = 0.5;
    int sign1;
        vector nnn[];
      foreach() {
            if (ff5[] <= 0. || ff5[] >= 1.) {
            foreach_dimension()
                    nnn.x[] = 0.;
            }
            else {
            coord m = mycs (point, ff5);
            foreach_dimension()
                    nnn.x[] = m.x;
            }
      }
    for(int phase=0;phase<=1;phase++){
   //     printf("pid=%d, T_leon 36, phase=%d\n",pid(),phase);
	sign1=2*phase-1;
        // foreach(noauto){
            foreach(){
           // if( ((phase==1 && ff5[]>=threshold) || (phase==0 && (1.0-ff5[])>threshold)) && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1))){
	 //  if(level==level_interface && css_test3_n[]>0.0 && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1))){
        //  if(level==level_interface && css_test3_n[]>=1.0 && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1)) && (!is_boundary(cell))){ //20230603
        if(level==level_interface && css_test3_n[]>=0.5 && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1)) && (!is_boundary(cell))){ //20231106
                double colinearity=0.0;
                coord pure;
                pure.x=x,pure.y=y; //,pure.z=z;
                double puretoint=HUGE;
                foreach_neighbor(2){ //at least 6 layers near the interface
                  //  if(level==level_interface && topo_mask[]==0 && css_test3_n[]>0.0 && (ff5[]>1e-6) && ((1.0-ff5[])>1e-6) && (!is_boundary(cell))){
                    // if(level==level_interface && topo_mask[]==0 && css_test3_n[]>=1.0 && (ff5[]>1e-6) && ((1.0-ff5[])>1e-6) && (!is_boundary(cell))){
                      if(level==level_interface && topo_mask[]==0 && css_test3_n[]>=0.5 && (ff5[]>1e-6) && ((1.0-ff5[])>1e-6) && (!is_boundary(cell))){
                        //calculate colinearity
                            //origin on interface cell
                            double colinearity_temp;
                            // coord nn = interface_normal7(point,ff5,hh2);
                            // coord nn =  mycs (point, ff5); //2023050603
                            coord nn;
                            foreach_dimension(){
                                nn.x = nnn.x[];
                            }
                            double alpha2=plane_alpha(ff5[],nn);;
                            coord rela_areacp,real_areacp;
                            plane_area_center (nn, alpha2,&rela_areacp);
			    normalize(&nn);
                            real_areacp.x = x + Delta*rela_areacp.x;
                            real_areacp.y = y + Delta*rela_areacp.y;
                           // real_areacp.z = z + Delta*rela_areacp.z;
                            coord mixed;
                            coord deltapure;
                            double leng1;
                            mixed.x=x,mixed.y=y;//,mixed.z=z;
                            foreach_dimension(){
                                deltapure.x = pure.x - mixed.x;
                            } 
                            leng1 = sqrt(sq(deltapure.x)+sq(deltapure.y));//+sq(deltapure.z));
                            if(leng1<0.0001){
                              //  printf("pid=%d leng1=%g\n",pid(),leng1);
                            }
                            //colinearity_temp = fabs(nn.x*deltapure.x/leng1+nn.y*deltapure.y/leng1+nn.z*deltapure.z/leng1);
                            colinearity_temp = fabs(nn.x*deltapure.x/leng1+nn.y*deltapure.y/leng1);//+nn.z*deltapure.z/leng1);

                            if(colinearity_temp>colinearity){
                                colinearity = colinearity_temp;
                                //calculate distance from interface to pure cell
                                    //position of interface center = real_areacp
                                    //distance between interface center and pure cell;
                                    puretoint = fabs(nn.x*(pure.x-real_areacp.x)+nn.y*(pure.y-real_areacp.y));//+nn.z*(pure.z-real_areacp.z));
                                   // distanceT[]=puretoint;
                                    
                            }
                    }
                }

                if(puretoint < (Delta*3.0*sqrt(2.0))){ //20230103 add
	                if(puretoint > Delta*0.1){
                            distanceT[] = puretoint;
                            double temp8 = (T[]-Tsat00)/distanceT[];
                            // if(fabs(temp8)>1000.0 && new_flag){
                            //    // printf("pid=%d,temp8:%g,%g,%g,%g,%g,%g\n",pid(),temp8,ff5[],css_test3[],css_test3_n[],T[],distanceT[]);
                            // }
                            if(phase==0){
                                phaseg[]=temp8;
                            }else{
                                phasel[]=-temp8;
                                // phasel[]=temp8; //20230918
                            }
                  }else{
                            distanceT[] = Delta*0.1;
                            double temp8 = (T[]-Tsat00)/distanceT[];
                            // if(fabs(temp8)>1000.0 && new_flag){
                            //    // printf("pid=%d,temp8:%g,%g,%g,%g,%g,%g\n",pid(),temp8,ff5[],css_test3[],css_test3_n[],T[],distanceT[]);
                            // }
                            if(phase==0){
                                phaseg[]=temp8;
                            }else{
                                phasel[]=-temp8;
                                // phasel[]=temp8; //20230918
                            }
		                             printf("puretoint less than 0.1*Delta, puretoint=%g ff=%g globali=%d\n", puretoint,ff5[],globali);
		             }
                }

            }

        }
    }
     //   printf("pid=%d, T_leon halfway\n",pid());
   //calulate Tgradient for interface cells
  for(int phase=0;phase<=1;phase++){
    sign1 = 2*phase-1;
    // foreach(noauto){
        foreach(){
           get_novalue_flag[]=0;
          //  if(level==level_interface && topo_mask[]==0 && css_test3_n[]>0.0 &&  (ff5[]>1e-6) && ((1.0-ff5[])>1e-6) ){
           if(level==level_interface && topo_mask[]==0 && css_test3_n[]>0.0 &&  (ff5[]>1e-6) && ((1.0-ff5[])>1e-6) &&  (!is_boundary(cell))){
                coord mixed;
                mixed.x=x,mixed.y=y;//,mixed.z=z;
                // coord nn=interface_normal7(point,ff5,hh2);;
                 coord nn =  mycs (point, ff5); //2023050603
                normalize(&nn);

                double weight_tot=0.0;
                double val_tot=0.0;
                int num_tot=0;

                double gradT_inner=0.0;
                foreach_neighbor(2){
                   // if( ((phase==1 && ff5[]>=threshold) || (phase==0 && (1.0-ff5[])>threshold)) && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1))){
		  // if(level==level_interface && css_test3_n[]>0.0 && (!is_boundary(cell)) && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1))){
          //  if(level==level_interface && css_test3_n[]>=1.0 && (!is_boundary(cell)) && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1)) && (!is_boundary(cell))){
            if(level==level_interface && css_test3_n[]>=0.5 && (!is_boundary(cell)) && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1)) && (!is_boundary(cell))){
                        coord pure;
                        pure.x = x, pure.y = y;//, pure.z = z;
                        coord puretomixed;
                        puretomixed.x = pure.x - mixed.x;
                        puretomixed.y = pure.y - mixed.y;
                       // puretomixed.z = pure.z - mixed.z;
                        double leng1 = sqrt(sq(puretomixed.x)+sq(puretomixed.y));//+sq(puretomixed.z));
                        double colinearity;
                       // if(fabs(leng1)>0.0001){
                               // colinearity = fabs(nn.x*puretomixed.x/leng1+nn.y*puretomixed.y/leng1+nn.z*puretomixed.z/leng1);
                          //3d      colinearity = fabs(nn.x*puretomixed.x/leng1+nn.y*puretomixed.y/leng1+nn.z*puretomixed.z/leng1)*(1.0/(leng1/Delta+1e-30));
                                // colinearity = fabs(nn.x*puretomixed.x/leng1+nn.y*puretomixed.y/leng1);//*(1.0/(leng1/Delta+1e-30));

                                // colinearity = fabs(nn.x*puretomixed.x/leng1+nn.y*puretomixed.y/leng1)*(1.0/(leng1/Delta+1e-30)); // before 20230614


                                colinearity = fabs(nn.x*puretomixed.x/leng1+nn.y*puretomixed.y/leng1);
                                //colinearity = fabs(nn.x*puretomixed.x/leng1+nn.y*puretomixed.y/leng1)*(1.0/(leng1/Delta+1e-30));
                                if(distanceT[]<(Delta*3.0*sqrt(2.0))){
                    //   if(colinearity < ){
                                    weight_tot += colinearity;
                                    if(phase==1){
                                        val_tot += phasel[]*colinearity;
                                    }else{
                                        val_tot += phaseg[]*colinearity;
                                    }
                                    num_tot += 1;
                                    // if(fabs(x-0.507812)<0.001 && fabs(y-0.382812)<0.001 && fabs(z-0.304688)<0.001 && new_flag){
                                    //     printf("phaseg[]=%g,phasel[]=%g,colinearity=%g\n",phaseg[],phasel[],colinearity);
                                    // }
                //	  }
                                }
                      //  }
                        
                    }
                } //foreach_neighbor
                if(num_tot>=1){
                        if(phase==1){
                             if(weight_tot<0.001){
                               // printf("weight_tot=%g val_tot=%g phasel=%g\n",weight_tot,val_tot,phasel[]);
                            }
                            phasel[]=val_tot/weight_tot;
                            // if(fabs(phasel[])>1000.0 && new_flag ){
                            //    // printf("inside leon: phasel=%g,weight=%g,ff=%g\n",phasel[],weight_tot,ff5[]);
                            //     foreach_neighbor(){
                            //        // printf("     T=%g,%g,%g,%g,%g\n",T[],x,y,z,ff5[]);
                            //     }
                            // }
			   
                        }else{
                            if(weight_tot<0.001){
                               //printf("weight_tot=%g val_tot=%g phaseg=%g\n",weight_tot,val_tot,phaseg[]);
                            }
                            phaseg[]=val_tot/weight_tot;
			                // if(fabs(phaseg[])>1000.0 && new_flag){
                            //    // printf("inside leon:phaseTgrad=%g,weight=%g,ff=%g,%g,%g,%g\n",phaseg[],weight_tot,ff5[],x,y,z);
                            //     foreach_neighbor(){
                            //      //   printf("     T=%g,%g,%g,%g,%g\n",T[],x,y,z,ff5[]);
                            //     }
                            // }
                        }
                }else{ //for microlayer, when the fluid interface is attached to solid cell, in between no liquid cell
                    get_novalue_flag[]=1;
                    coord m;
                    bool flag = false;

                    bool intersect_true_neighbour=false;
                    foreach_neighbor(1){
                        if(intersect_true[]==1){
                            intersect_true_neighbour = true;
                        }
                    }

                    // if(intersect_true[]==1 || intersect_true_neighbour){
                    //  if(1==1 && intersect_true[]!=1){
                    if(1==1){
                            if(phase==0){
                                // if(css_test2[]>css_test[]){
                                if(ff5[]<=0.5){
                                    // m = mycs (point, css_test2);
                                    // double alpha_g=plane_alpha(css_test2[],m);
                                    m = mycs (point, ff_oppo);
                                    double alpha_g=plane_alpha(ff_oppo[],m);
                                    coord area_center;
                                    plane_area_center (m, alpha_g, &area_center);
                                    double length_g = sqrt(sq(area_center.x) + sq(area_center.y));
                                     if(length_g>length_threshold){
                                        phaseg[] = (T[] - Tsat00)/(length_g*Delta);
                                        flag = true;
                                    }else{
                                         phaseg[] = (T[] - Tsat00)/(length_threshold*Delta);
                                        flag = true;
                                    }
                                }else{  //get average T[] from neighbors
                                    if(1==1){
                                        // m = mycs (point, css_test2);
                                        // double alpha_g=plane_alpha(css_test2[],m);
                                        m = mycs (point, ff_oppo);
                                        double alpha_g=plane_alpha(ff_oppo[],m);
                                        coord area_center;
                                        plane_area_center (m, alpha_g, &area_center);
                                        double length_g = sqrt(sq(area_center.x) + sq(area_center.y));


                                        double total=0.0;
                                        double weight=0.0;
                                        foreach_neighbor(2){
                                            if(css_test3_n[]>0.5 && level==level_interface && ff5[]<=0.5){
                                                    total += T[];
                                                    weight +=1;
                                            }
                                        }
                                        if(weight>=1){
                                            double T_average = total/weight;
                                            if(length_g>length_threshold){
                                                phaseg[] = (T_average - Tsat00)/(length_g*Delta);
                                                flag = true;
                                            }else{
                                                phaseg[] = (T_average - Tsat00)/(length_threshold*Delta);
                                                flag = true;
                                            }

                                        }else{
                                            phaseg[] = 0.0;
                                        }

                                    }
                                        
                                        // flag = true;
                                }
                                if((!flag) && (aimg[]>=Tsat00)){
                                    if(css_test3[]<1.0 && css_test3[]>0.0){
                                         // m = mycs (point, css_test2);
                                        // double alpha_g=plane_alpha(css_test2[],m);
                                        m = mycs (point, ff_oppo);
                                         double alpha_g=plane_alpha(ff_oppo[],m);
                                        coord area_center;
                                        plane_area_center (m, alpha_g, &area_center);
                                        
                                        coord ms = mycs (point, css_test3);
                                        double alpha_s=plane_alpha(css_test3[],ms);
                                        coord area_centers;
                                        plane_area_center (ms, alpha_s, &area_centers);
                                        double length_sg = sqrt(sq(area_center.x-area_centers.x) + sq(area_center.y-area_centers.y));
                                         if(length_sg>length_threshold){
                                            phaseg[] = (aimg[] - Tsat00)/(Delta*length_sg);
                                        }else{
                                            phaseg[] = (aimg[] - Tsat00)/(Delta*length_threshold);
                                        }
                                    }
                                }
                            }else{ //phase==1
                                // if(css_test[]>=css_test2[]){
                                if(ff5[]>=0.5){
                                    // m = mycs (point, css_test);
                                    // double alpha_l=plane_alpha(css_test[],m);
                                    m = mycs (point, ff5);
                                    double alpha_l=plane_alpha(ff5[],m);
                                    coord area_center;
                                    plane_area_center (m, alpha_l, &area_center);
                                    double length_l = sqrt(sq(area_center.x) + sq(area_center.y));
                                    if(length_l>length_threshold){
                                        phasel[] = -(T[] - Tsat00)/(length_l*Delta);
                                        flag = true;
                                    }else{
                                         phasel[] = -(T[] - Tsat00)/(length_threshold*Delta);
                                        flag = true;
                                    }
                                }else{//get average T[] from neighbors
                                    

                                    if(1==1){
                                        // m = mycs (point, css_test);
                                        // double alpha_l=plane_alpha(css_test[],m);
                                        m = mycs (point, ff5);
                                        double alpha_l=plane_alpha(ff5[],m);
                                        coord area_center;
                                        plane_area_center (m, alpha_l, &area_center);
                                        double length_l = sqrt(sq(area_center.x) + sq(area_center.y));


                                        double total=0.0;
                                        double weight=0.0;
                                        foreach_neighbor(2){
                                            if(css_test3_n[]>0.5 && level==level_interface && ff5[]>=0.5){
                                                    total += T[];
                                                    weight +=1;
                                            }
                                        }
                                        if(weight>=1){
                                            double T_average = total/weight;
                                            if(length_l>length_threshold){
                                                phasel[] = -(T_average - Tsat00)/(length_l*Delta);
                                                flag = true;
                                            }else{
                                                phasel[] = -(T_average - Tsat00)/(length_threshold*Delta);
                                                flag = true;
                                            }

                                        }else{
                                            phasel[] = 0.0;
                                        }

                                    }
                                        // flag = true;
                                }
                                if((!flag) && (aiml[]>=Tsat00)){
                                    if(css_test3[]<1.0 && css_test3[]>0.0){
                                       // m = mycs (point, css_test);
                                        // double alpha_l=plane_alpha(css_test[],m);
                                        m = mycs (point, ff5);
                                        double alpha_l=plane_alpha(ff5[],m);
                                        coord area_center;
                                        plane_area_center (m, alpha_l, &area_center);
                                        
                                        coord ms = mycs (point, css_test3);
                                        double alpha_s=plane_alpha(css_test3[],ms);
                                        coord area_centers;
                                        plane_area_center (ms, alpha_s, &area_centers);
                                        double length_sl = sqrt(sq(area_center.x-area_centers.x) + sq(area_center.y-area_centers.y));
                                        if(length_sl>length_threshold){
                                            phasel[] = -(aiml[] - Tsat00)/(Delta*length_sl);
                                        }else{
                                            phasel[] = -(aiml[] - Tsat00)/(Delta*length_threshold);
                                        }
                                    }
                                }
                            }
                    }
                }
            }//topo_mask
    }
    if(((q.aiml.i && phase==1) ||  (q.aimg.i && phase==0)) && 1==0){
        // printf("dalaji hahaha\n");
      vector nnn_s[];
      vector center_s[];
      foreach() {
                coord m;
                if (css_test3[] <= 0. || css_test3[] >= 1.) {
                  foreach_dimension(){
                          m.x = 0;
                          nnn_s.x[] = 0.;
                  }
                }
                else {
                  m = mycs (point, css_test3);
                  foreach_dimension()
                          nnn_s.x[] = m.x;
                  double alpha2=plane_alpha(css_test3[],m);;
                  coord center_s1;
                  plane_area_center (m, alpha2,&center_s1);
                  foreach_dimension(){
                      center_s.x[] = center_s1.x;
                  }
                }
      }
        foreach(){
          if(fabs(get_novalue_flag[]-1)<1e-6){
              // printf("get_novalue_flag[]==1 aiml[]=%g\n",q.aiml[]);
              coord nn; //direction of liquid
              foreach_dimension(){
                            nn.x = nnn.x[];
              }
              double alpha2=plane_alpha(ff5[],nn);;
              coord center_f;
              plane_area_center (nn, alpha2,&center_f);
              
              // if(pid()==0){
              //     printf("q.aiml.i exist and enter\n");
              //         fprintf(stderr,"q.aiml.i exist and enter\n");
              // }
              double distance_select=HUGE;
              double value=0.0;
              foreach_neighbor(1){ //could include itself, triple point
                  if(css_test3[]>0.0 && css_test3[]<1.0){
                      if(phase==1){
                            double direction_juge = nn.x*nnn_s.x[] +nn.y*nnn_s.y[];
                            if(fabs(q.aiml[])>0.0 && direction_juge>0.0){
                                    double distance = sqrt(sq(center_f.x-center_s.x[])+sq(center_f.y-center_s.y[]));
                                    // distance = max(0.1,distance)*Delta;
                                    if(distance_select>distance){
                                      distance_select=distance;
                                      value = (Tsat00 - q.aiml[]); //-temp8;
                                    }
                              }
                      }else{
                            double direction_juge = -nn.x*nnn_s.x[] +(-nn.y)*nnn_s.y[];
                            if(fabs(q.aimg[])>0.0 && direction_juge>0.0){
                                    double distance = sqrt(sq(center_f.x-center_s.x[])+sq(center_f.y-center_s.y[]));
                                    // distance = max(0.1,distance)*Delta;
                                    if(distance_select>distance){
                                      distance_select=distance;
                                      value = (q.aimg[] - Tsat00); //temp8;
                                    }
                              }
                      }
                  }
              }
              if(distance_select<HUGE/2.0){
                    // distance_select = max(0.1,distance_select)*Delta;
                    distance_select = max(length_threshold,distance_select)*Delta;
                    if(phase==1){
                        phasel[] = value/distance_select;
                    }else{
                        phaseg[] = value/distance_select;
                    }
              }
              //  printf("get_novalue_flag[]==1 aiml[]=%g phasel[]=%g distance_select=%g\n",q.aiml[],phasel[],distance_select);
              
          }
        }
    }


  } //phase
   //   printf("pid=%d, T_leon finish\n",pid());
   if (is_constant_m) {
      foreach(){
        phasel[] = -38.5;
        phaseg[] = 0.0;
      }
   }
} 



// void Tgrad_leon(scalar ff5,scalar phaseg,scalar phasel){
// void Tgrad_leon_combined(struct Tgradleon q){
void Tgrad_leon(struct Tgradleon q){
  scalar ff5 = q.ff5;
  scalar phaseg = q.phaseg;
  scalar phasel = q.phasel;
  double length_threshold = 0.05;//0.01; 
 bool is_constant_m=false;//;true;//false;//true;
    //ff5 should be original ff
 //       printf("pid=%d, T_leon begin\n",pid());
    vector hh2[];
 //   printf("pid=%d, T_leon 22\n",pid());
    heights(ff5,hh2);
 //   printf("pid=%d, T_leon 24\n",pid());
    scalar distanceT[];
    scalar get_novalue_flag[];
  //1: calculate Tgradient for pure cells
      //require topo_mask
    foreach(){
        phaseg[]=0.0;
        phasel[]=0.0;
        distanceT[] = HUGE;
    }
    double threshold = 0.5;
    int sign1;
        vector nnn[];
      foreach() {
            if (ff5[] <= 0. || ff5[] >= 1.) {
            foreach_dimension()
                    nnn.x[] = 0.;
            }
            else {
            coord m = mycs (point, ff5);
            foreach_dimension()
                    nnn.x[] = m.x;
            }
      }
      foreach(){
        ff_oppo[] = 1.0-ff5[];
      }
    for(int phase=0;phase<=1;phase++){
   //     printf("pid=%d, T_leon 36, phase=%d\n",pid(),phase);
	sign1=2*phase-1;
        // foreach(noauto){
            foreach(){
           // if( ((phase==1 && ff5[]>=threshold) || (phase==0 && (1.0-ff5[])>threshold)) && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1))){
	 //  if(level==level_interface && css_test3_n[]>0.0 && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1))){
        //  if(level==level_interface && css_test3_n[]>=1.0 && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1)) && (!is_boundary(cell))){ //20230603
        if(level==level_interface && css_test3_n[]>=0.5 && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1)) && (!is_boundary(cell))){ //20231106
                double colinearity=0.0;
                coord pure;
                pure.x=x,pure.y=y; //,pure.z=z;
                double puretoint=HUGE;
                foreach_neighbor(2){ //at least 6 layers near the interface
                  //  if(level==level_interface && topo_mask[]==0 && css_test3_n[]>0.0 && (ff5[]>1e-6) && ((1.0-ff5[])>1e-6) && (!is_boundary(cell))){
                    // if(level==level_interface && topo_mask[]==0 && css_test3_n[]>=1.0 && (ff5[]>1e-6) && ((1.0-ff5[])>1e-6) && (!is_boundary(cell))){
                      if(level==level_interface && topo_mask[]==0 && css_test3_n[]>=0.5 && (ff5[]>1e-6) && ((1.0-ff5[])>1e-6) && (!is_boundary(cell))){
                        //calculate colinearity
                            //origin on interface cell
                            double colinearity_temp;
                            // coord nn = interface_normal7(point,ff5,hh2);
                            // coord nn =  mycs (point, ff5); //2023050603
                            coord nn;
                            foreach_dimension(){
                                nn.x = nnn.x[];
                            }
                            double alpha2=plane_alpha(ff5[],nn);;
                            coord rela_areacp,real_areacp;
                            plane_area_center (nn, alpha2,&rela_areacp);
			    normalize(&nn);
                            real_areacp.x = x + Delta*rela_areacp.x;
                            real_areacp.y = y + Delta*rela_areacp.y;
                           // real_areacp.z = z + Delta*rela_areacp.z;
                            coord mixed;
                            coord deltapure;
                            double leng1;
                            mixed.x=x,mixed.y=y;//,mixed.z=z;
                            foreach_dimension(){
                                deltapure.x = pure.x - mixed.x;
                            } 
                            leng1 = sqrt(sq(deltapure.x)+sq(deltapure.y));//+sq(deltapure.z));
                            if(leng1<0.0001){
                              //  printf("pid=%d leng1=%g\n",pid(),leng1);
                            }
                            //colinearity_temp = fabs(nn.x*deltapure.x/leng1+nn.y*deltapure.y/leng1+nn.z*deltapure.z/leng1);
                            colinearity_temp = fabs(nn.x*deltapure.x/leng1+nn.y*deltapure.y/leng1);//+nn.z*deltapure.z/leng1);

                            if(colinearity_temp>colinearity){
                                colinearity = colinearity_temp;
                                //calculate distance from interface to pure cell
                                    //position of interface center = real_areacp
                                    //distance between interface center and pure cell;
                                    puretoint = fabs(nn.x*(pure.x-real_areacp.x)+nn.y*(pure.y-real_areacp.y));//+nn.z*(pure.z-real_areacp.z));
                                   // distanceT[]=puretoint;
                                    
                            }
                    }
                }

                if(puretoint < (Delta*3.0*sqrt(2.0))){ //20230103 add
	                if(puretoint > Delta*length_threshold){
                            distanceT[] = puretoint;
                            
                            // if(fabs(temp8)>1000.0 && new_flag){
                            //    // printf("pid=%d,temp8:%g,%g,%g,%g,%g,%g\n",pid(),temp8,ff5[],css_test3[],css_test3_n[],T[],distanceT[]);
                            // }
                            if(phase==0){
                                double aa = Tg[];
                                double temp8 = (aa-Tsat00)/distanceT[];
                                phaseg[]=temp8;
                            }else{
                                double aa = Tl[];
                                double temp8 = (aa-Tsat00)/distanceT[];
                                phasel[]=-temp8;
                                // phasel[]=temp8; //20230918
                            }
                  }else{
                            distanceT[] = Delta*length_threshold;
                            
                            // if(fabs(temp8)>1000.0 && new_flag){
                            //    // printf("pid=%d,temp8:%g,%g,%g,%g,%g,%g\n",pid(),temp8,ff5[],css_test3[],css_test3_n[],T[],distanceT[]);
                            // }
                            if(phase==0){
                                double aa = Tg[];
                                double temp8 = (aa-Tsat00)/distanceT[];
                                phaseg[]=temp8;
                            }else{
                                double aa = Tl[];
                                double temp8 = (aa-Tsat00)/distanceT[];
                                phasel[]=-temp8;
                                // phasel[]=temp8; //20230918
                            }
		                             printf("puretoint less than 0.01*Delta, puretoint=%g ff=%g globali=%d\n", puretoint,ff5[],globali);
		             }
                }

            }

        }
    }
     //   printf("pid=%d, T_leon halfway\n",pid());
   //calulate Tgradient for interface cells
  for(int phase=0;phase<=1;phase++){
    sign1 = 2*phase-1;
    // foreach(noauto){
        foreach(){
           get_novalue_flag[]=0;
          //  if(level==level_interface && topo_mask[]==0 && css_test3_n[]>0.0 &&  (ff5[]>1e-6) && ((1.0-ff5[])>1e-6) ){
           if(level==level_interface && topo_mask[]==0 && css_test3_n[]>0.0 &&  (ff5[]>1e-6) && ((1.0-ff5[])>1e-6) &&  (!is_boundary(cell))){
                coord mixed;
                mixed.x=x,mixed.y=y;//,mixed.z=z;
                // coord nn=interface_normal7(point,ff5,hh2);;
                 coord nn =  mycs (point, ff5); //2023050603
                normalize(&nn);

                double weight_tot=0.0;
                double val_tot=0.0;
                int num_tot=0;

                double gradT_inner=0.0;

                int num_neighbors=0;
                 double gradients[25], distances[25];

                foreach_neighbor(2){
                   // if( ((phase==1 && ff5[]>=threshold) || (phase==0 && (1.0-ff5[])>threshold)) && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1))){
		  // if(level==level_interface && css_test3_n[]>0.0 && (!is_boundary(cell)) && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1))){
          //  if(level==level_interface && css_test3_n[]>=1.0 && (!is_boundary(cell)) && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1)) && (!is_boundary(cell))){
            if(level==level_interface && css_test3_n[]>=0.5 && (!is_boundary(cell)) && ((topo_mask[]==1*sign1) || (topo_mask[]==2*sign1)) && (!is_boundary(cell))){
                        coord pure;
                        pure.x = x, pure.y = y;//, pure.z = z;
                        coord puretomixed;
                        puretomixed.x = pure.x - mixed.x;
                        puretomixed.y = pure.y - mixed.y;
                       // puretomixed.z = pure.z - mixed.z;
                        double leng1 = sqrt(sq(puretomixed.x)+sq(puretomixed.y));//+sq(puretomixed.z));
                        double colinearity;
                       // if(fabs(leng1)>0.0001){
                               // colinearity = fabs(nn.x*puretomixed.x/leng1+nn.y*puretomixed.y/leng1+nn.z*puretomixed.z/leng1);
                          //3d      colinearity = fabs(nn.x*puretomixed.x/leng1+nn.y*puretomixed.y/leng1+nn.z*puretomixed.z/leng1)*(1.0/(leng1/Delta+1e-30));
                                // colinearity = fabs(nn.x*puretomixed.x/leng1+nn.y*puretomixed.y/leng1);//*(1.0/(leng1/Delta+1e-30));

                                // colinearity = fabs(nn.x*puretomixed.x/leng1+nn.y*puretomixed.y/leng1)*(1.0/(leng1/Delta+1e-30)); // before 20230614


                                colinearity = fabs(nn.x*puretomixed.x/leng1+nn.y*puretomixed.y/leng1);
                                //colinearity = fabs(nn.x*puretomixed.x/leng1+nn.y*puretomixed.y/leng1)*(1.0/(leng1/Delta+1e-30));
                                if(distanceT[]<(Delta*3.0*sqrt(2.0))){
                    //   if(colinearity < ){
                                    weight_tot += colinearity;
                                    if(phase==1){
                                        val_tot += phasel[]*colinearity;
                                        gradients[num_neighbors] = phasel[]; 
                                        distances[num_neighbors] = distanceT[]/2.0;
                                        num_neighbors += 1;
                                    }else{
                                        val_tot += phaseg[]*colinearity;
                                        gradients[num_neighbors] = phaseg[]; 
                                        distances[num_neighbors] = distanceT[]/2.0;
                                        num_neighbors += 1;
                                    }
                                    num_tot += 1;
                                    
                                    
                                    // if(fabs(x-0.507812)<0.001 && fabs(y-0.382812)<0.001 && fabs(z-0.304688)<0.001 && new_flag){
                                    //     printf("phaseg[]=%g,phasel[]=%g,colinearity=%g\n",phaseg[],phasel[],colinearity);
                                    // }
                //	  }
                                }
                      //  }
                        
                    }
                } //foreach_neighbor
                if(num_tot>=1){
                    bool flag=false;
                    // if(num_neighbors >= 3){
                    //     bool tooCloseForQuadratic = are_points_too_close(distances, num_neighbors, 0.5*Delta);
                    //     if(!tooCloseForQuadratic){
                    //         double a, b_quad, c; // Parameters for quadratic regression
                    //         double interface_gradient_quadratic;
                    //         if (quadratic_regression(distances, gradients, num_neighbors, &a, &b_quad, &c)) {
                    //             interface_gradient_quadratic = predict_quadratic(a, b_quad, c, 0.0);
                    //             flag = true;
                    //             // printf("Estimated temperature gradient at the interface (Quadratic): %f\n", interface_gradient_quadratic);
                    //         } else {
                    //             // printf("Quadratic regression failed.\n");
                    //         }
                    //         if(flag){
                    //             if(phase==1){
                    //                 phasel[]=interface_gradient_quadratic;
                    //             }else{
                    //                 phaseg[]=interface_gradient_quadratic;
                    //             }
                    //         }
                    //     }
                    // }
                    // if(!flag && num_neighbors>=2){
                    //     bool tooCloseForQuadratic = are_points_too_close(distances, num_neighbors, 0.5*Delta);
                    //     if(!tooCloseForQuadratic){
                    //         double m, b; // Parameters for linear regression
                    //         double interface_gradient_linear;
                    //         if(linear_regression(distances, gradients, num_neighbors, &m, &b)){
                    //             interface_gradient_linear = predict_linear(m, b, 0.0);
                    //             flag = true;
                    //         }
                    //         if(flag){
                    //             if(phase==1){
                    //                 phasel[]=interface_gradient_linear;
                    //             }else{
                    //                 phaseg[]=interface_gradient_linear;
                    //             }
                    //         }
                    //     }
                    // }
                    if(!flag){
                        if(phase==1){
                             if(weight_tot<0.001){
                               // printf("weight_tot=%g val_tot=%g phasel=%g\n",weight_tot,val_tot,phasel[]);
                            }
                            phasel[]=val_tot/weight_tot;
                            // if(fabs(phasel[])>1000.0 && new_flag ){
                            //    // printf("inside leon: phasel=%g,weight=%g,ff=%g\n",phasel[],weight_tot,ff5[]);
                            //     foreach_neighbor(){
                            //        // printf("     T=%g,%g,%g,%g,%g\n",T[],x,y,z,ff5[]);
                            //     }
                            // }
			   
                        }else{
                            if(weight_tot<0.001){
                               //printf("weight_tot=%g val_tot=%g phaseg=%g\n",weight_tot,val_tot,phaseg[]);
                            }
                            phaseg[]=val_tot/weight_tot;
			                // if(fabs(phaseg[])>1000.0 && new_flag){
                            //    // printf("inside leon:phaseTgrad=%g,weight=%g,ff=%g,%g,%g,%g\n",phaseg[],weight_tot,ff5[],x,y,z);
                            //     foreach_neighbor(){
                            //      //   printf("     T=%g,%g,%g,%g,%g\n",T[],x,y,z,ff5[]);
                            //     }
                            // }
                        }
                    }
                }else{ //for microlayer, when the fluid interface is attached to solid cell, in between no liquid cell
                    get_novalue_flag[]=1;
                    coord m;
                    bool flag = false;

                    bool intersect_true_neighbour=false;
                    foreach_neighbor(1){
                        if(intersect_true[]==1){
                            intersect_true_neighbour = true;
                        }
                    }

                    // if(intersect_true[]==1 || intersect_true_neighbour){
                    //  if(1==1 && intersect_true[]!=1){
                    if(1==1){
                            if(phase==0){
                                // if(css_test2[]>css_test[]){
                                // if(ff5[]<=0.5){
                                if(ff5[]<1.0){
                                    // m = mycs (point, css_test2);
                                    // double alpha_g=plane_alpha(css_test2[],m);
                                    // coord area_center;
                                    // plane_area_center (m, alpha_g, &area_center);
                                    // double length_g = sqrt(sq(area_center.x) + sq(area_center.y));
                                    double length_g = css_test2[]/2.0;
                                     if(length_g>length_threshold){
                                        phaseg[] = (Tg[] - Tsat00)/(length_g*Delta);
                                        flag = true;
                                    }else{
                                         phaseg[] = (Tg[] - Tsat00)/(length_threshold*Delta);
                                        flag = true;
                                    }
                                } // }else{  //get average T[] from neighbors
                                //     if(1==1){
                                //         m = mycs (point, css_test2);
                                //         double alpha_g=plane_alpha(css_test2[],m);
                                //         coord area_center;
                                //         plane_area_center (m, alpha_g, &area_center);
                                //         double length_g = sqrt(sq(area_center.x) + sq(area_center.y));


                                //         double total=0.0;
                                //         double weight=0.0;
                                //         foreach_neighbor(2){
                                //             if(css_test3_n[]>0.5 && level==level_interface && ff5[]<=0.5){
                                //                     total += T[];
                                //                     weight +=1;
                                //             }
                                //         }
                                //         if(weight>=1){
                                //             double T_average = total/weight;
                                //             if(length_g>length_threshold){
                                //                 phaseg[] = (T_average - Tsat00)/(length_g*Delta);
                                //                 flag = true;
                                //             }else{
                                //                 phaseg[] = (T_average - Tsat00)/(length_threshold*Delta);
                                //                 flag = true;
                                //             }

                                //         }else{
                                //             phaseg[] = 0.0;
                                //         }

                                //     }
                                        
                                //         // flag = true;
                                // }
                                if((!flag) && (aimg[]>=Tsat00)){
                                    if(css_test3[]<1.0 && css_test3[]>0.0){
                                        // m = mycs (point, css_test2);
                                        // double alpha_g=plane_alpha(css_test2[],m);
                                        m = mycs (point, ff_oppo);
                                        double alpha_g=plane_alpha(ff_oppo[],m);
                                        coord area_center;
                                        plane_area_center (m, alpha_g, &area_center);
                                        
                                        coord ms = mycs (point, css_test3);
                                        double alpha_s=plane_alpha(css_test3[],ms);
                                        coord area_centers;
                                        plane_area_center (ms, alpha_s, &area_centers);
                                        double length_sg = sqrt(sq(area_center.x-area_centers.x) + sq(area_center.y-area_centers.y));
                                         if(length_sg>length_threshold){
                                            phaseg[] = (aimg[] - Tsat00)/(Delta*length_sg);
                                        }else{
                                            phaseg[] = (aimg[] - Tsat00)/(Delta*length_threshold);
                                        }
                                    }
                                }
                            }else{ //phase==1
                                // if(css_test[]>=css_test2[]){
                                // if(ff5[]>=0.5){
                                 if(ff5[]>0){
                                    // m = mycs (point, css_test);
                                    // double alpha_l=plane_alpha(css_test[],m);
                                    // coord area_center;
                                    // plane_area_center (m, alpha_l, &area_center);
                                    // double length_l = sqrt(sq(area_center.x) + sq(area_center.y));
                                    double length_l = css_test[]/2.0;
                                    if(length_l>length_threshold){
                                        phasel[] = -(Tl[] - Tsat00)/(length_l*Delta);
                                        flag = true;
                                    }else{
                                         phasel[] = -(Tl[] - Tsat00)/(length_threshold*Delta);
                                        flag = true;
                                    }
                                }
                                // }else{//get average T[] from neighbors
                                    

                                //     if(1==1){
                                //         m = mycs (point, css_test);
                                //         double alpha_l=plane_alpha(css_test[],m);
                                //         coord area_center;
                                //         plane_area_center (m, alpha_l, &area_center);
                                //         double length_l = sqrt(sq(area_center.x) + sq(area_center.y));


                                //         double total=0.0;
                                //         double weight=0.0;
                                //         foreach_neighbor(2){
                                //             if(css_test3_n[]>0.5 && level==level_interface && ff5[]>=0.5){
                                //                     total += T[];
                                //                     weight +=1;
                                //             }
                                //         }
                                //         if(weight>=1){
                                //             double T_average = total/weight;
                                //             if(length_l>length_threshold){
                                //                 phasel[] = -(T_average - Tsat00)/(length_l*Delta);
                                //                 flag = true;
                                //             }else{
                                //                 phasel[] = -(T_average - Tsat00)/(length_threshold*Delta);
                                //                 flag = true;
                                //             }

                                //         }else{
                                //             phasel[] = 0.0;
                                //         }

                                //     }
                                //         // flag = true;
                                // }
                                if((!flag) && (aiml[]>=Tsat00)){
                                    if(css_test3[]<1.0 && css_test3[]>0.0){
                                        // m = mycs (point, css_test);
                                        // double alpha_l=plane_alpha(css_test[],m);
                                        m = mycs (point, ff5);
                                        double alpha_l=plane_alpha(ff5[],m);
                                        coord area_center;
                                        plane_area_center (m, alpha_l, &area_center);
                                        
                                        coord ms = mycs (point, css_test3);
                                        double alpha_s=plane_alpha(css_test3[],ms);
                                        coord area_centers;
                                        plane_area_center (ms, alpha_s, &area_centers);
                                        double length_sl = sqrt(sq(area_center.x-area_centers.x) + sq(area_center.y-area_centers.y));
                                        if(length_sl>length_threshold){
                                            phasel[] = -(aiml[] - Tsat00)/(Delta*length_sl);
                                        }else{
                                            phasel[] = -(aiml[] - Tsat00)/(Delta*length_threshold);
                                        }
                                    }
                                }
                            }
                    }
                }
            }//topo_mask
    }
    if(((q.aiml.i && phase==1) ||  (q.aimg.i && phase==0)) && 1==0){
        // printf("dalaji hahaha\n");
      vector nnn_s[];
      vector center_s[];
      foreach() {
                coord m;
                if (css_test3[] <= 0. || css_test3[] >= 1.) {
                  foreach_dimension(){
                          m.x = 0;
                          nnn_s.x[] = 0.;
                  }
                }
                else {
                  m = mycs (point, css_test3);
                  foreach_dimension()
                          nnn_s.x[] = m.x;
                  double alpha2=plane_alpha(css_test3[],m);;
                  coord center_s1;
                  plane_area_center (m, alpha2,&center_s1);
                  foreach_dimension(){
                      center_s.x[] = center_s1.x;
                  }
                }
      }
        foreach(){
          if(fabs(get_novalue_flag[]-1)<1e-6){
              // printf("get_novalue_flag[]==1 aiml[]=%g\n",q.aiml[]);
              coord nn; //direction of liquid
              foreach_dimension(){
                            nn.x = nnn.x[];
              }
              double alpha2=plane_alpha(ff5[],nn);;
              coord center_f;
              plane_area_center (nn, alpha2,&center_f);
              
              // if(pid()==0){
              //     printf("q.aiml.i exist and enter\n");
              //         fprintf(stderr,"q.aiml.i exist and enter\n");
              // }
              double distance_select=HUGE;
              double value=0.0;
              foreach_neighbor(1){ //could include itself, triple point
                  if(css_test3[]>0.0 && css_test3[]<1.0){
                      if(phase==1){
                            double direction_juge = nn.x*nnn_s.x[] +nn.y*nnn_s.y[];
                            if(fabs(q.aiml[])>0.0 && direction_juge>0.0){
                                    double distance = sqrt(sq(center_f.x-center_s.x[])+sq(center_f.y-center_s.y[]));
                                    // distance = max(0.1,distance)*Delta;
                                    if(distance_select>distance){
                                      distance_select=distance;
                                      value = (Tsat00 - q.aiml[]); //-temp8;
                                    }
                              }
                      }else{
                            double direction_juge = -nn.x*nnn_s.x[] +(-nn.y)*nnn_s.y[];
                            if(fabs(q.aimg[])>0.0 && direction_juge>0.0){
                                    double distance = sqrt(sq(center_f.x-center_s.x[])+sq(center_f.y-center_s.y[]));
                                    // distance = max(0.1,distance)*Delta;
                                    if(distance_select>distance){
                                      distance_select=distance;
                                      value = (q.aimg[] - Tsat00); //temp8;
                                    }
                              }
                      }
                  }
              }
              if(distance_select<HUGE/2.0){
                    distance_select = max(length_threshold,distance_select)*Delta;
                    if(phase==1){
                        phasel[] = value/distance_select;
                    }else{
                        phaseg[] = value/distance_select;
                    }
              }
              //  printf("get_novalue_flag[]==1 aiml[]=%g phasel[]=%g distance_select=%g\n",q.aiml[],phasel[],distance_select);
              
          }
        }
    }


  } //phase
//   scalar grad_embed_3rd_g[],grad_embed_3rd_l[];
  //gas
  foreach(){
    ff_oppo[] = 1.0-ff5[];
    T_grad_method_g[]=0;
    T_grad_method_l[]=0;
    grad_embed_3rd_g[] = 0;
    grad_embed_3rd_l[] = 0;
  }
  foreach(){
    if(css_test3[]<1.0 && level==level_interface){
        if(ff_oppo[]<1.0-1e-6 && ff_oppo[]>1e-6){
            grad_embed_3rd_g[] = 0;
            double method=0;
            double grad=0.0;
            double bc = Tsat00;
            grad = embed_3_order_Tgrad_gas (point, Tg, ff_oppo, &method, bc);
            T_grad_method_g[]=method;
            if(method>=1){
                grad_embed_3rd_g[] = -grad;
            }
        }
    }
  }
  foreach(){
    if(css_test3[]<1.0 && level==level_interface){
        if(ff5[]<1.0-1e-6 && ff5[]>1e-6){
            grad_embed_3rd_l[] = 0;
            double method=0;
            double grad=0.0;
            double bc = Tsat00;
            grad = embed_3_order_Tgrad_liquid (point, Tl, ff5, &method, bc);
            T_grad_method_l[]=method;
            if(method>=1){
                grad_embed_3rd_l[] = grad;
            }
        }
    }
  }
  if(1==0){
        foreach(){
            if(css_test3[]<1.0 && level==level_interface){
                if(ff5[]<1.0-1e-6 && ff5[]>1e-6){
                    if(T_grad_method_g[]>=2){ // no first order
                        phaseg[]=grad_embed_3rd_g[];
                    }else{
                        T_grad_method_g[]=4;
                    }
                    if(T_grad_method_l[]>=2){ //
                        phasel[]=grad_embed_3rd_l[];
                    }else{
                        T_grad_method_l[]=4;
                    }
                }
            }
        }
  }
} 

//input css_test3_n solid; the direction of phase0and1 is the normal from (gas and water) to solid
// but with horizontal solid we don't need this
// void Tgrad_leon_s(scalar ff5,scalar phaseg,scalar phasel){
void Tgrad_leon_s(struct Tgradleon q){
  scalar ff5 = q.ff5;
  scalar phaseg = q.phaseg;
  scalar phasel = q.phasel;
    //ff5 should be original ff
        // printf("pid=%d, T_leon begin\n",pid());
    vector hh2[];
    // printf("pid=%d, T_leon 22\n",pid());
    heights(ff5,hh2);
    // printf("pid=%d, T_leon 24\n",pid());
    scalar distanceT[];
  //1: calculate Tgradient for pure cells
      //require topo_mask_s
    foreach(){
        phaseg[]=0.0;
        phasel[]=0.0;
        distanceT[] = HUGE;
    }
    double threshold = 0.5;
    int sign1;
    for(int phase=0;phase<=1;phase++){
    //    printf("pid=%d, T_leon 36, phase=%d\n",pid(),phase);
	sign1=2*phase-1;
        foreach(noauto){
           // if( ((phase==1 && ff5[]>=threshold) || (phase==0 && (1.0-ff5[])>threshold)) && ((topo_mask_s[]==1*sign1) || (topo_mask_s[]==2*sign1))){
	   if(level==level_interface && ((topo_mask_s[]==1*sign1) || (topo_mask_s[]==2*sign1))){
                double colinearity=0.0;
                coord pure;
                pure.x=x,pure.y=y;//,pure.z=z;
                double puretoint=HUGE;
                foreach_neighbor(2){ //at least 6 layers near the interface
                    if(level==level_interface && topo_mask_s[]==0 && (ff5[]>1e-6) && ((1.0-ff5[])>1e-6) && (!is_boundary(cell))){
                        //calculate colinearity
                            //origin on interface cell
                            double colinearity_temp;
                            coord nn = interface_normal7(point,ff5,hh2);
                            double alpha2=plane_alpha(ff5[],nn);;
                            coord rela_areacp,real_areacp;
                            plane_area_center (nn, alpha2,&rela_areacp);
			    normalize(&nn);
                            real_areacp.x = x + Delta*rela_areacp.x;
                            real_areacp.y = y + Delta*rela_areacp.y;
                          //  real_areacp.z = z + Delta*rela_areacp.z;
                            coord mixed;
                            coord deltapure;
                            double leng1;
                            mixed.x=x,mixed.y=y,mixed.z=z;
                            foreach_dimension(){
                                deltapure.x = pure.x - mixed.x;
                            } 
                            leng1 = sqrt(sq(deltapure.x)+sq(deltapure.y));//+sq(deltapure.z));
                            if(leng1<0.0001){
                              //  printf("pid=%d leng1=%g\n",pid(),leng1);
                            }
                            colinearity_temp = fabs(nn.x*deltapure.x/leng1+nn.y*deltapure.y/leng1);//+nn.z*deltapure.z/leng1);

                            if(colinearity_temp>colinearity){
                                colinearity = colinearity_temp;
                                //calculate distance from interface to pure cell
                                    //position of interface center = real_areacp
                                    //distance between interface center and pure cell;
                                    puretoint = fabs(nn.x*(pure.x-real_areacp.x)+nn.y*(pure.y-real_areacp.y));//+nn.z*(pure.z-real_areacp.z));
                                   // distanceT[]=puretoint;
                                    
                            }
                    }
                }

                if(puretoint < (Delta*3.0*sqrt(2.0))){ //20230103 add
	                if(puretoint > Delta*0.1){
		                distanceT[] = puretoint;
                    double temp8 = (T[]-Tsat00)/distanceT[];
                    if(phase==0){
                        phaseg[]=temp8;
                    }else{
                        phasel[]=-temp8;
                    }
                  }else{
		                             printf("puretoint less than 0.1*Delta, puretoint=%g ff=%g globali=%d\n", puretoint,ff5[],globali);
                      // distanceT[] = Delta*0.1;//puretoint;
                      // double temp8 = (T[]-Tsat00)/distanceT[];
                      // if(phase==0){
                      //   phaseg[]=temp8;
                      // }else{
                      //     phasel[]=-temp8;
                      // }
		             }
                }

            }

        }
    }
        printf("pid=%d, T_leon halfway\n",pid());
   //calulate Tgradient for interface cells
  for(int phase=0;phase<=1;phase++){
    sign1 = 2*phase-1;
    foreach(noauto){
            if(level==level_interface && topo_mask_s[]==0 &&  (ff5[]>1e-6) && ((1.0-ff5[])>1e-6)){
                coord mixed;
                mixed.x=x,mixed.y=y;//,mixed.z=z;
                coord nn=interface_normal7(point,ff5,hh2);;
                normalize(&nn);

                double weight_tot=0.0;
                double val_tot=0.0;
                int num_tot=0;
                foreach_neighbor(2){
                   // if( ((phase==1 && ff5[]>=threshold) || (phase==0 && (1.0-ff5[])>threshold)) && ((topo_mask_s[]==1*sign1) || (topo_mask_s[]==2*sign1))){
		   if(level==level_interface && (!is_boundary(cell)) && ((topo_mask_s[]==1*sign1) || (topo_mask_s[]==2*sign1))){
                        coord pure;
                        pure.x = x, pure.y = y;//, pure.z = z;
                        coord puretomixed;
                        puretomixed.x = pure.x - mixed.x;
                        puretomixed.y = pure.y - mixed.y;
                   //     puretomixed.z = pure.z - mixed.z;
                        double leng1 = sqrt(sq(puretomixed.x)+sq(puretomixed.y));//+sq(puretomixed.z));
                        double colinearity;
                        colinearity = fabs(nn.x*puretomixed.x/leng1+nn.y*puretomixed.y/leng1);//+nn.z*puretomixed.z/leng1);
                        if(distanceT[]<(Delta*3.0*sqrt(2.0))){
		       //   if(colinearity < ){
                            weight_tot += colinearity;
                            if(phase==1){
                                val_tot += phasel[]*colinearity;
                            }else{
                                val_tot += phaseg[]*colinearity;
                            }
                            num_tot += 1;
		//	  }
                        }
                        
                    }
                } //foreach_neighbor
                if(num_tot>=1){
                        if(phase==1){
                             if(weight_tot<0.001){
                                printf("weight_tot=%g val_tot=%g phasel=%g\n",weight_tot,val_tot,phasel[]);
                            }
                            phasel[]=val_tot/weight_tot;
			   
                        }else{
                            if(weight_tot<0.001){
                              printf("weight_tot=%g val_tot=%g phaseg=%g\n",weight_tot,val_tot,phaseg[]);
                            }
                            phaseg[]=val_tot/weight_tot;
			    
                        }
                }
            }//topo_mask_s
    }
  } //phase
      printf("pid=%d, T_leon finish\n",pid());
}




bool linear_regression(double *x, double *y, int n, double *m, double *b) {
    double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
    for (int i = 0; i < n; i++) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xy += x[i] * y[i];
        sum_x2 += x[i] * x[i];
    }

    double denominator = (n * sum_x2 - sum_x * sum_x);
    if (denominator == 0) {
        return false; // Prevent division by zero
    }

    *m = (n * sum_xy - sum_x * sum_y) / denominator;
    *b = (sum_y - *m * sum_x) / n;
    return true;
}

bool quadratic_regression(double *x, double *y, int n, double *a, double *b, double *c) {
    if (n < 3) return false;

    double sum_x = 0, sum_x2 = 0, sum_x3 = 0, sum_x4 = 0, sum_y = 0, sum_xy = 0, sum_x2y = 0;
    for (int i = 0; i < n; i++) {
        double xi = x[i], yi = y[i];
        sum_x += xi;
        sum_x2 += xi * xi;
        sum_x3 += xi * xi * xi;
        sum_x4 += xi * xi * xi * xi;
        sum_y += yi;
        sum_xy += xi * yi;
        sum_x2y += xi * xi * yi;
    }

    double D = n * (sum_x2 * sum_x4 - sum_x3 * sum_x3) - sum_x * (sum_x * sum_x4 - sum_x3 * sum_x2) + sum_x2 * (sum_x * sum_x3 - sum_x2 * sum_x2);
    if (D == 0) return false;

    *a = (sum_y * (sum_x2 * sum_x4 - sum_x3 * sum_x3) - sum_x * (sum_xy * sum_x4 - sum_x3 * sum_x2y) + sum_x2 * (sum_xy * sum_x3 - sum_x2 * sum_x2y)) / D;
    *b = (n * (sum_xy * sum_x4 - sum_x3 * sum_x2y) - sum_y * (sum_x * sum_x4 - sum_x2 * sum_x3) + sum_x2 * (sum_x * sum_x2y - sum_xy * sum_x2)) / D;
    *c = (n * (sum_x2 * sum_x2y - sum_xy * sum_x3) - sum_x * (sum_x * sum_x2y - sum_x2 * sum_xy) + sum_y * (sum_x * sum_x3 - sum_x2 * sum_x2)) / D;

    return true;
}

double predict_linear(double m, double b, double x) {
    return m * x + b;
}

double predict_quadratic(double a, double b, double c, double x) {
    return a * x * x + b * x + c;
}

bool are_points_too_close(double *distances, int n, double delta) {
    double min_distance = distances[0];
    double max_distance = distances[0];
    for (int i = 1; i < n; i++) {
        if (distances[i] < min_distance) min_distance = distances[i];
        if (distances[i] > max_distance) max_distance = distances[i];
    }
    return (max_distance - min_distance) <= delta;
}
