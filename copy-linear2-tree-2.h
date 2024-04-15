#include "fractions.h"

// distance from pp to interface


//double point_to_line(coord nn, double alpha_n, coord pp){
//     double val1=fabs(nn.x*pp.x+nn.y*pp.y-alpha_n)/fabs(nn.x*nn.x+nn.y*nn.y);
//     return val1;
//}

//#define CCCCMIN 0.01

// artribute should be refined before scalar s
double lim_cut = 1e-12;
double point_to_line(coord nn, double alpha_n, coord pp){
    if(sqrt(nn.x*nn.x+nn.y*nn.y+nn.z*nn.z)<lim_cut){
         return sqrt(3.0)/2.0;
     }

     double a = nn.x/sqrt(nn.x*nn.x+nn.y*nn.y+nn.z*nn.z), b=nn.y/sqrt(nn.x*nn.x+nn.y*nn.y+nn.z*nn.z);
     double c = nn.z/sqrt(nn.x*nn.x+nn.y*nn.y+nn.z*nn.z);
     double alpha_nn = alpha_n/sqrt(nn.x*nn.x+nn.y*nn.y+nn.z*nn.z);
     double val1=fabs(a*pp.x+b*pp.y+c*pp.z-alpha_nn)/fabs(sqrt(a*a+b*b+c*c));
    
    if(val1<=sqrt(3.0)/2.0)
      return val1;
    else if(val1>sqrt(3.0)/2.0){
      return sqrt(3.0)/2.0;
    }
     //return 1.0;
}


#include "curvature.h"
coord interface_normal3 (Point point, scalar c, vector h)
{
  coord n;
  if (!h.x.i || (n = height_normal (point, c, h)).x == nodata)
    n = mycs (point, c);
  return n;
}


attribute {
   scalar ff3; //volume of fluid
   //face vector ffs2;
   vector centroid_p; //center of the cell. 
   //vector height_temp;
  // vector normal;  //this normal using interface_nomral3 to calculate normal;
   double Tsat0;
   double T_wall;

   scalar centroidpx;
   scalar centroidpy;
   scalar centroidpz;

   vector hhh;

   
}
extern double L0;
extern int isboundary_right;
extern double boundary_value_right;


void Tl_refine2(Point point, scalar s){  //for T_l

    

     scalar ff3 = s.ff3;



     double cc = ff3[];
     double ss = s[];
     double T_wall = s.T_wall;
     vector centroid_p = s.centroid_p;

     bool flag_wallT = false;
     if(isboundary_right==1 && is_boundary(neighbor(1,0))){
            flag_wallT = true;
     }
     
   if(cc <= 0 ){
     foreach_child(){
        s[] = Tsat0*0.0; 
     }
   }else{
     
     
     double Delta_coarse = Delta;

       //coord n = mycs(point,ff3);  // this normal could be changed to being caculate by height later
       coord n = interface_normal3(point,ff3,hhh); //normal got from height function.
       double alpha = plane_alpha (cc, n);

      //for alpha_refine
        coord m;
        double alphac = 2.*alpha;
        foreach_dimension()
              m.x = n.x;

     foreach_child(){
      //  if((coarse(ff3) > 0. && coarse(ff3,child.x) > 0. &&
      //			   coarse(ff3,0,child.y) > 0. && coarse(ff3,child.x,child.y) > 0.)){
            
      //   }
         // 2-dimensions
        coord p_c_child_c; // tansfer to parent's coordate

        static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
            coord nc;
            foreach_dimension(){
                 nc.x = child.x * n.x;
            }
            double cc2;
            cc2 = rectangle_fraction(nc,alpha,a,b);  //child_ff3
            bool flag_child_0 = false;
           double alpha_child=alphac;
      if(cc2>0.0 && cc2<1.0){  //if centroid_p refined before s, then there is no need to recaculate p_c_child_c
            //double alpha_child=alphac;
            foreach_dimension()
                    alpha_child -= child.x*m.x/2.;  //child_alpha
        
            coord child_p;
            coord child_normal = n;
            plane_center(child_normal,alpha_child,cc2,&child_p);
           // foreach_dimension(){
                     p_c_child_c.x = child_p.x/2.0 + child.x*0.25;
                     p_c_child_c.y = child_p.y/2.0 + child.y*0.25;
           // }
       }else if(cc2>=1){
           p_c_child_c.x = child.x*0.25;
           p_c_child_c.y = child.y*0.25;
       }else{ //cc2<=0
           s[] = Tsat0*0.0;
           flag_child_0 = true;
       }

  if(!flag_child_0){
     double weight_i;
     int ilim = 2*child.x;
     int jlim = 2*child.y;

     double weight_tol=0.0;
     int weight_n=0;
     double val_tol=0.0;
     for(int i=0;abs(i)<abs(ilim);i+=child.x){
        for(int j=0;abs(j)<abs(jlim);j+=child.y){
           bool flag_i = true;
               if(coarse(centroid_p.x,i,j)>HUGE/2.0){
                   flag_i = false ;
                }
               if(coarse(centroid_p.y,i,j)>HUGE/2.0){
                   flag_i = false ;
               }
         if(flag_i){
               
               double dd = sqrt(sq(p_c_child_c.x-coarse(centroid_p.x,i,j)) + sq(p_c_child_c.y-coarse(centroid_p.y,i,j)));
               if (dd>1e-12){
                      weight_i = 1.0/dd;
               }else{
                      weight_i = 1.0/(1e-12);
               }
               weight_tol += weight_i;
               weight_n++;
               val_tol += weight_i*coarse(s,i,j);
         }
        }
      }//for child.

      // interface_affect
     /*
      if(cc2>0.0 && cc2<1.0){
          // double dd = fabs(point_to_line(child_normal,alpha_child,p_c_child_c)); //non-dimensional
          double dd = fabs(point_to_line(n,alpha_child,p_c_child_c)); //non-dimensional
          dd = dd/2.0;
          if (dd>1e-12){
                      weight_i = 1.0/dd;
          }else{
                      weight_i = 1.0/(1e-12);
          }
          weight_tol += weight_i;
          weight_n++;
          val_tol += weight_i*coarse(s);
      }
     */
     /*
     if(flag_wallT){ // right_temperature_wall
         double dd = fabs(L0-coarse(centroid_p.x));
         dd = dd/2.0;
         if (dd>1e-12){
                      weight_i = 1.0/dd;
          }else{
                      weight_i = 1.0/(1e-12);
          }
          weight_tol += weight_i;
          weight_n++;
          val_tol += weight_i*T_wall;
     }
     */

         if(weight_n>=1){
            s[] = val_tol/weight_tol*ff3[];
         }else{
            s[] = Tsat0*ff3[];
         }
       } // !!flag_child_0

          
     }//foreach_child


   }

}




void Tl_restriction2(Point point, scalar s){
  scalar ff3 = s.ff3;
  //double cc = ff3[];
  double sum = 0.;
 // int part_cell=0;
  foreach_child(){
    if(ff3[]>0){
     //part_cell++;
     sum += ff3[]*s[];  //since 
    }
  }
  s[] = sum/(1 << dimension)/(ff3[] + 1e-30)*ff3[];
}






//the following version is for both inverse and !inverse,
// but 
// need height, restriction(ff),poisson_phase,Tsat0,Twall,inverse,restriction(centroid_p,centroidpx,centroidpy);
void Tl_refine33(Point point, scalar s){  //for T_l

    // write for ff3 need to change when add leafs, so ff3 should be ff only;
    //beside I should write an refine and restriction for centroid_p_l and centroid_p_g
     double Tsat0 = s.Tsat0;
     scalar ff3 = s.ff3;

     double cc;

     double weight_limit =1e-12;



  //   char name153[80];
  //   sprintf(name153,"Tl_refine33_childff_check%g.dat",t);
  //   FILE * fp153 = fopen(name153,"a");
  //   fprintf(fp153,"parent: %g %g\n", ff3[],EPS);
  //  // if(ff3[]>EPS && ff3[]<1-EPS){
  //        foreach_child(){
  //            if(ff3[]>EPS && ff3[]<1-EPS)
  //               fprintf(fp153,"%g %g %g !!!!!!!!!!!!!!!!!!!\n",x,y,ff3[]);
  //        }
  //  // }
  //   fclose(fp153);



if(s.inverse){
     cc = 1.0-ff3[];  
}else{
     cc = ff3[];
}

     double ss = s[];
     double T_wall = s.T_wall;
     vector centroid_p = s.centroid_p;

     scalar centroidpx = s.centroidpx;
     scalar centroidpy = s.centroidpy;
     scalar centroidpz = s.centroidpz;

     bool flag_wallT = false;
     if(isboundary_right==1 && is_boundary(neighbor(1,0))){
            flag_wallT = true;
     }
     
   if(cc <= 0 ){
     foreach_child(){
        s[] = Tsat0; 
     }
   }else{
     
    
     double Delta_coarse = Delta;

       //coord n = mycs(point,ff3);  // this normal could be changed to being caculate by height later
       coord n = interface_normal3(point,ff3,hhh); //normal got from height function.
if(s.inverse){
    foreach_dimension(){
       n.x = -n.x;
    }
}
       double alpha = plane_alpha (cc, n);

      //for alpha_refine
        coord m;
        double alphac = 2.*alpha;
        foreach_dimension()
              m.x = n.x;

     foreach_child(){
      //  if((coarse(ff3) > 0. && coarse(ff3,child.x) > 0. &&
      //			   coarse(ff3,0,child.y) > 0. && coarse(ff3,child.x,child.y) > 0.)){
            
      //   }
         // 2-dimensions
        
        coord p_c_child_c; // tansfer to parent's coordate
        double cc2;
        double alpha_child;
        bool flag_child_0 = false;

     if(cc<1){    
  
        static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
            coord nc;
            foreach_dimension(){
                 nc.x = child.x * n.x;
            }
            
            cc2 = rectangle_fraction(nc,alpha,a,b);  //child_ff3
            
           alpha_child=alphac;
          if(cc2>0.0 && cc2<1.0){  //if centroid_p refined before s, then there is no need to recaculate p_c_child_c
                //double alpha_child=alphac;
                foreach_dimension()
                        alpha_child -= child.x*m.x/2.;  //child_alpha
            
                coord child_p;
                coord child_normal = n;
                plane_center(child_normal,alpha_child,cc2,&child_p);
              // foreach_dimension(){
                        p_c_child_c.x = child_p.x/2.0 + child.x*0.25;
                        p_c_child_c.y = child_p.y/2.0 + child.y*0.25;
                        p_c_child_c.z = child_p.z/2.0 + child.z*0.25;
              // }
          }else if(cc2>=1){
              p_c_child_c.x = child.x*0.25;
              p_c_child_c.y = child.y*0.25;
              p_c_child_c.z = child.z*0.25;
          }else{ //cc2<=0
              s[] = Tsat0;
              flag_child_0 = true;
          }
     }else{
          cc2 = 1.0;
           p_c_child_c.x = child.x*0.25;
           p_c_child_c.y = child.y*0.25;
           p_c_child_c.z = child.z*0.25;
     }
  if(!flag_child_0){
     double weight_i;
     int ilim = 2*child.x;
     int jlim = 2*child.y;
     int klim = 2*child.z;

      double weight_tol=0.0;
     int weight_n=0;
     double val_tol=0.0;
     bool too_close1=false;
     for(int i=0;abs(i)<abs(ilim);i+=child.x){
        for(int j=0;abs(j)<abs(jlim);j+=child.y){
            for(int k=0;abs(k)<abs(klim);k+=child.z){
              if(!too_close1){         
                        // bool flag_i = true;
                        bool skip1=false;
                        double ssss;
                            if(!s.inverse){
                                  /*
                        if(coarse(centroid_p.x,i,j)>HUGE/2.0){
                            flag_i = false ;
                          }
                        if(coarse(centroid_p.y,i,j)>HUGE/2.0){
                            flag_i = false ;
                        }
                                  */
                                ssss= coarse(ff3,i,j,k);  
                            }else{
                                /*
                                    if(fabs(coarse(centroid_p.x,i,j)-0.0)<1e-20){
                            flag_i = false ;
                          }
                            if(!s.inverse){  
                                   p_temp.x = coarse(centroidpx,i,j);
                                   p_temp.y = coarse(centroidpy,i,j);
                        }      if(fabs(coarse(centroid_p.y,i,j)-0.0)<1e-20){
                            flag_i = false ;
                        }
                                */
                                ssss = 1.0 - coarse(ff3,i,j,k); 
                            }
                            if(ssss<=EPS){
                                skip1 = true;        
                            }
                            
                      if(!skip1){
                          coord p_temp;
                          if(!s.inverse){  
                            // double dd = sqrt(sq(p_c_child_c.x-coarse(centroid_p.x,i,j)) + sq(p_c_child_c.y-coarse(centroid_p.y,i,j)));
                                  //p_temp.x = coarse(centroid_p.x,i,j);
                                  //p_temp.y = coarse(centroid_p.y,i,j);
                                  p_temp.x = coarse(centroidpx,i,j,k);
                                  p_temp.y = coarse(centroidpy,i,j,k);
                                  p_temp.z = coarse(centroidpz,i,j,k);                                  
                          }else{
                            //if(coarse(centroid_p.x,i,j)>HUGE/2.0){
                            if(ssss>=1.0){
                                  p_temp.x = 0.0;
                                  p_temp.y = 0.0;
                                  p_temp.z = 0.0;
                            }else{
                        double sss = coarse(ff3,i,j,k);
                                    if(1.0-sss < EPS){
                                            // skip1 = true;
                                            p_temp.x = 0.0;
                                            p_temp.y = 0.0;
                                            p_temp.z = 0.0;
                                    }else{
                                            //p_temp.x = -sss*coarse(centroid_p.x,i,j)/(1.0-sss);
                          //p_temp.y = -sss*coarse(centroid_p.y,i,j)/(1.0-sss);
                                            p_temp.x = -sss*coarse(centroidpx,i,j,k)/(1.0-sss);
                                            p_temp.y = -sss*coarse(centroidpy,i,j,k)/(1.0-sss);
                                            p_temp.z = -sss*coarse(centroidpz,i,j,k)/(1.0-sss);
                                            foreach_dimension(){
                                                  if(p_temp.x>=0.5-EPS || p_temp.x<=-0.5+EPS){
                                                    skip1 = true;
                                                  }
                                            }
                                    }
                       
                            }
                          }

                            if(!skip1){
                                p_temp.x = i + p_temp.x;
                                p_temp.y = j + p_temp.y;
                                p_temp.z = k + p_temp.z;
                                //double dd = sqrt(sq(p_c_child_c.x-coarse(centroid_p.x,i,j)) + sq(p_c_child_c.y-coarse(centroid_p.y,i,j)));
                                double dd = sqrt(sq(p_c_child_c.x-p_temp.x) + sq(p_c_child_c.y-p_temp.y) + sq(p_c_child_c.z-p_temp.z));
                                if (dd>weight_limit){
                                        weight_i = 1.0/dd;
                                      if(ssss>EPS){
                                            weight_tol += weight_i;
                                            weight_n++;
                                            val_tol += weight_i*coarse(s,i,j,k);  //20221003

                                                  // char name154[80];
                                                  // sprintf(name154,"ill_point_refine33_check%g.dat",t);
                                                  // FILE * fp154 = fopen(name154,"a");
                                                  // fprintf(fp154,"first: %g %g %g %g %g\n",x,y,weight_i,ssss,coarse(s,i,j)/ssss);
                                                  // fclose(fp154); 
                                      }
                                }else{
                                        
                                        if(ssss>EPS){
                                                too_close1 = true;
                                                s[] = coarse(s,i,j,k);

                                                  // char name154[80];
                                                  // sprintf(name154,"ill_point_refine33_check%g.dat",t);
                                                  // FILE * fp154 = fopen(name154,"a");
                                                  // fprintf(fp154,"final2: %g %g %g %g\n",x,y,weight_i,coarse(s,i,j)/ssss);
                                                  // fclose(fp154); 
                                        }
                                        //else{
                                         //    s[] = Tsat0*cc2;
                                        //}
                                        //weight_i = 1.0/(weight_limit);
                                }
                               
                            }// compelete1
                      }
              }//too_close
          } //k
        }  // j
      }//for child. i 

      // interface_affect
 

//  char name154[80];
//  sprintf(name154,"ill_point_refine33_check%g.dat",t);
//  FILE * fp154 = fopen(name154,"a");
//  fprintf(fp154,"weight_n: %g %g %d\n",x,y,weight_n);
//  fclose(fp154); 



bool interface_effect = true;

         // if(!too_close1){
                bool too_close_to_interface = false;
if(interface_effect){
                //if(cc2>0.0 && cc2<1.0){
                if(cc2>=EPS && (cc2<=1-EPS) ){
                     coord temp_p = (coord){0.0,0.0,0.0};
                    // double dd = fabs(point_to_line(n,alpha_child,p_c_child_c)); //non-dimensional
                     double dd = fabs(point_to_line(n,alpha_child,temp_p)); //non-dimensional
                    //double dd = fabs(point_to_line(n,alpha,p_c_child_c)); //non-dimensional
                   dd = dd/2.0; // transfer distance to parent coodinate
                    if (dd>weight_limit){
                                weight_i = 1.0/dd;
                    }else{
                                //weight_i = 1.0/(weight_limit);
                                too_close_to_interface = true;
                                
                    }
                    if(!too_close_to_interface){
                        weight_tol += weight_i;
                        weight_n++;
                      // val_tol += weight_i*coarse(s);
                        val_tol += weight_i*Tsat0;
                    }
                }
    } // interface_effect           
                if((!too_close_to_interface)){
                   if(!too_close1){
                      if(flag_wallT){ // right_temperature_wall
                          if(child.y==-1){
                          double dd = fabs(0.5-coarse(centroid_p.x));
                          // dd = dd/2.0; // transfer distance to parent coodinate
                          if (dd>weight_limit){
                                        weight_i = 1.0/dd;
                            }else{
                                        weight_i = 1.0/(weight_limit);
                            }
                            weight_tol += weight_i;
                            weight_n++;
                            val_tol += weight_i*T_wall;
                          }
                          
                      }
                  

                      if(weight_n>=1){
                        // s[] = val_tol/weight_tol*ff3[];
                          // if(weight_tol>weight_limit)
                                  s[] = val_tol/weight_tol;
                          // else
                          //       s[] = Tsat0*cc2;
                                        //if(x>3.7 && x<3.74 && y>3.7 && y<3.74 ){
                                                  // char name154[80];
                                                  // sprintf(name154,"ill_point_refine33_check%g.dat",t);
                                                  // FILE * fp154 = fopen(name154,"a");
                                                  // fprintf(fp154,"final: %g %g %g %g %g\n",x,y,weight_tol,val_tol,s[]);
                                                  // fclose(fp154);    
                                       // }

                            if(s[]<0){
                                printf("weight_tol=%g\n",weight_tol);
                               // if(x>3.7 && x<3.74 && y>3.7 && y<3.74 ){
                                                  // char name154[80];
                                                  // sprintf(name154,"ill_point_refine33_check%g.dat",t);
                                                  // FILE * fp154 = fopen(name154,"a");
                                                  // fprintf(fp154,"suprising, s<0 !!??\n");
                                                  // fclose(fp154);    
                              //  }

                            }
                      }else{
                        //  s[] = Tsat0*ff3[];
                          s[] = Tsat0;
                        // s[] = Tsat0*0;

                                                  // char name154[80];
                                                  // sprintf(name154,"ill_point_refine33_check%g.dat",t);
                                                  // FILE * fp154 = fopen(name154,"a");
                                                  // fprintf(fp154,"final3: %g %g %g\n",s[],cc2,Tsat0);
                                                  // fclose(fp154);
                      }
                   }

                }else{
                     s[] =  Tsat0;
                                                  // char name154[80];
                                                  // sprintf(name154,"ill_point_refine33_check%g.dat",t);
                                                  // FILE * fp154 = fopen(name154,"a");
                                                  // fprintf(fp154,"final4: %g %g %g\n",s[],cc2,Tsat0);
                                                  // fclose(fp154);
                }


          // }//!too_close
       } // !!flag_child_0
       


     }//foreach_child


   }

 
}

void Tl_refine33_adap(Point point, scalar s){  //for T_l  // checked for 3D

    // write for ff3 need to change when add leafs, so ff3 should be ff only;
    //beside I should write an refine and restriction for centroid_p_l and centroid_p_g
     double Tsat0 = s.Tsat0;
     scalar ff3 = s.ff3;

     double cc;

     double weight_limit =1e-12;



  //   char name153[80];
  //   sprintf(name153,"Tl_refine33_childff_check%g.dat",t);
  //   FILE * fp153 = fopen(name153,"a");
  //   fprintf(fp153,"parent: %g %g\n", ff3[],EPS);
  //  // if(ff3[]>EPS && ff3[]<1-EPS){
  //        foreach_child(){
  //            if(ff3[]>EPS && ff3[]<1-EPS)
  //               fprintf(fp153,"%g %g %g !!!!!!!!!!!!!!!!!!!\n",x,y,ff3[]);
  //        }
  //  // }
  //   fclose(fp153);



if(s.inverse){
     cc = 1.0-ff3[];  
}else{
     cc = ff3[];
}

     double ss = s[];
     double T_wall = s.T_wall;
     vector centroid_p = s.centroid_p;

     scalar centroidpx = s.centroidpx;
     scalar centroidpy = s.centroidpy;
     scalar centroidpz = s.centroidpz;

     bool flag_wallT = false;
     if(isboundary_right==1 && is_boundary(neighbor(1,0))){
            flag_wallT = true;
     }
     
   if(cc <= 0 ){
     foreach_child(){
        s[] = Tsat0; 
     }
   }else{
     
    
     double Delta_coarse = Delta;

       //coord n = mycs(point,ff3);  // this normal could be changed to being caculate by height later
       coord n = interface_normal3(point,ff3,hhh); //normal got from height function.
if(s.inverse){
    foreach_dimension(){
       n.x = -n.x;
    }
}
       double alpha = plane_alpha (cc, n);

      //for alpha_refine
        coord m;
        double alphac = 2.*alpha;
        foreach_dimension()
              m.x = n.x;

     foreach_child(){
      //  if((coarse(ff3) > 0. && coarse(ff3,child.x) > 0. &&
      //			   coarse(ff3,0,child.y) > 0. && coarse(ff3,child.x,child.y) > 0.)){
            
      //   }
         // 2-dimensions
        
        coord p_c_child_c; // tansfer to parent's coordate
        double cc2;
        double alpha_child;
        bool flag_child_0 = false;

     if(cc<1){    
  
        static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
            coord nc;
            foreach_dimension(){
                 nc.x = child.x * n.x;
            }
            
            cc2 = rectangle_fraction(nc,alpha,a,b);  //child_ff3
            
           alpha_child=alphac;
          if(cc2>0.0 && cc2<1.0){  //if centroid_p refined before s, then there is no need to recaculate p_c_child_c
                //double alpha_child=alphac;
                foreach_dimension()
                        alpha_child -= child.x*m.x/2.;  //child_alpha
            
                coord child_p;
                coord child_normal = n;
                plane_center(child_normal,alpha_child,cc2,&child_p);
              // foreach_dimension(){
                        p_c_child_c.x = child_p.x/2.0 + child.x*0.25;
                        p_c_child_c.y = child_p.y/2.0 + child.y*0.25;
                        p_c_child_c.z = child_p.z/2.0 + child.z*0.25;
              // }
          }else if(cc2>=1){
              p_c_child_c.x = child.x*0.25;
              p_c_child_c.y = child.y*0.25;
              p_c_child_c.z = child.z*0.25;
          }else{ //cc2<=0
              s[] = Tsat0;
              flag_child_0 = true;
          }
     }else{
          cc2 = 1.0;
           p_c_child_c.x = child.x*0.25;
           p_c_child_c.y = child.y*0.25;
           p_c_child_c.z = child.z*0.25;
     }
  if(!flag_child_0){
     double weight_i;
     int ilim = 2*child.x;
     int jlim = 2*child.y;
     int klim = 2*child.z;

      double weight_tol=0.0;
     int weight_n=0;
     double val_tol=0.0;
     bool too_close1=false;
     for(int i=0;abs(i)<abs(ilim);i+=child.x){
        for(int j=0;abs(j)<abs(jlim);j+=child.y){
            for(int k=0;abs(k)<abs(klim);k+=child.z){
              if(!too_close1){         
                        // bool flag_i = true;
                        bool skip1=false;
                        double ssss;
                            if(!s.inverse){
                                  /*
                        if(coarse(centroid_p.x,i,j)>HUGE/2.0){
                            flag_i = false ;
                          }
                        if(coarse(centroid_p.y,i,j)>HUGE/2.0){
                            flag_i = false ;
                        }
                                  */
                                ssss= coarse(ff3,i,j,k);  
                            }else{
                                /*
                                    if(fabs(coarse(centroid_p.x,i,j)-0.0)<1e-20){
                            flag_i = false ;
                          }
                            if(!s.inverse){  
                                   p_temp.x = coarse(centroidpx,i,j);
                                   p_temp.y = coarse(centroidpy,i,j);
                        }      if(fabs(coarse(centroid_p.y,i,j)-0.0)<1e-20){
                            flag_i = false ;
                        }
                                */
                                ssss = 1.0 - coarse(ff3,i,j,k); 
                            }
                            if(ssss<=EPS){
                                skip1 = true;        
                            }
                            
                      if(!skip1){
                          coord p_temp;
                          if(!s.inverse){  
                            // double dd = sqrt(sq(p_c_child_c.x-coarse(centroid_p.x,i,j)) + sq(p_c_child_c.y-coarse(centroid_p.y,i,j)));
                                  //p_temp.x = coarse(centroid_p.x,i,j);
                                  //p_temp.y = coarse(centroid_p.y,i,j);
                                  p_temp.x = coarse(centroidpx,i,j,k);
                                  p_temp.y = coarse(centroidpy,i,j,k);
                                  p_temp.z = coarse(centroidpz,i,j,k);
                          }else{
                            //if(coarse(centroid_p.x,i,j)>HUGE/2.0){
                            if(ssss>=1.0){
                                  p_temp.x = 0.0;
                                  p_temp.y = 0.0;
                                  p_temp.z = 0.0;
                            }else{
                        double sss = coarse(ff3,i,j,k);
                                    if(1.0-sss < EPS){
                                            // skip1 = true;
                                            p_temp.x = 0.0;
                                            p_temp.y = 0.0;
                                            p_temp.z = 0.0;
                                    }else{
                                            //p_temp.x = -sss*coarse(centroid_p.x,i,j)/(1.0-sss);
                          //p_temp.y = -sss*coarse(centroid_p.y,i,j)/(1.0-sss);
                                            p_temp.x = -sss*coarse(centroidpx,i,j,k)/(1.0-sss);
                                            p_temp.y = -sss*coarse(centroidpy,i,j,k)/(1.0-sss);
                                            p_temp.z = -sss*coarse(centroidpz,i,j,k)/(1.0-sss);
                                            foreach_dimension(){
                                                  if(p_temp.x>=0.5-EPS || p_temp.x<=-0.5+EPS){
                                                    skip1 = true;
                                                  }
                                            }
                                    }
                       
                            }
                          }

                            if(!skip1){
                                p_temp.x = i + p_temp.x;
                                p_temp.y = j + p_temp.y;
                                p_temp.z = k + p_temp.z;
                                //double dd = sqrt(sq(p_c_child_c.x-coarse(centroid_p.x,i,j)) + sq(p_c_child_c.y-coarse(centroid_p.y,i,j)));
                                double dd = sqrt(sq(p_c_child_c.x-p_temp.x) + sq(p_c_child_c.y-p_temp.y) + sq(p_c_child_c.z-p_temp.z));
                                if (dd>weight_limit){
                                        weight_i = 1.0/dd;
                                      if(ssss>EPS){
                                            weight_tol += weight_i;
                                            weight_n++;
                                            val_tol += weight_i*coarse(s,i,j,k);  //20221003

                                                  // char name154[80];
                                                  // sprintf(name154,"ill_point_refine33_check%g.dat",t);
                                                  // FILE * fp154 = fopen(name154,"a");
                                                  // fprintf(fp154,"first: %g %g %g %g %g\n",x,y,weight_i,ssss,coarse(s,i,j)/ssss);
                                                  // fclose(fp154); 
                                      }
                                }else{
                                        
                                        if(ssss>EPS){
                                                too_close1 = true;
                                                s[] = coarse(s,i,j,k);

                                                  // char name154[80];
                                                  // sprintf(name154,"ill_point_refine33_check%g.dat",t);
                                                  // FILE * fp154 = fopen(name154,"a");
                                                  // fprintf(fp154,"final2: %g %g %g %g\n",x,y,weight_i,coarse(s,i,j)/ssss);
                                                  // fclose(fp154); 
                                        }
                                        //else{
                                         //    s[] = Tsat0*cc2;
                                        //}
                                        //weight_i = 1.0/(weight_limit);
                                }
                               
                            }// compelete1
                      }
              }//too_close
          }  //z
        }  //y
      }//for child.  //x

      // interface_affect
 

//  char name154[80];
//  sprintf(name154,"ill_point_refine33_check%g.dat",t);
//  FILE * fp154 = fopen(name154,"a");
//  fprintf(fp154,"weight_n: %g %g %d\n",x,y,weight_n);
//  fclose(fp154); 



bool interface_effect = false;

         // if(!too_close1){
                bool too_close_to_interface = false;
if(interface_effect){
                //if(cc2>0.0 && cc2<1.0){
                if(cc2>=EPS && (cc2<=1-EPS) ){
                     coord temp_p = (coord){0.0,0.0,0.0};
                    // double dd = fabs(point_to_line(n,alpha_child,p_c_child_c)); //non-dimensional
                     double dd = fabs(point_to_line(n,alpha_child,temp_p)); //non-dimensional
                    //double dd = fabs(point_to_line(n,alpha,p_c_child_c)); //non-dimensional
                   dd = dd/2.0; // transfer distance to parent coodinate
                    if (dd>weight_limit){
                                weight_i = 1.0/dd;
                    }else{
                                //weight_i = 1.0/(weight_limit);
                                too_close_to_interface = true;
                                
                    }
                    if(!too_close_to_interface){
                        weight_tol += weight_i;
                        weight_n++;
                      // val_tol += weight_i*coarse(s);
                        val_tol += weight_i*Tsat0;
                    }
                }
    } // interface_effect           
                if((!too_close_to_interface)){
                   if(!too_close1){
                      if(flag_wallT){ // right_temperature_wall
                          if(child.y==-1){
                          double dd = fabs(0.5-coarse(centroid_p.x));
                          // dd = dd/2.0; // transfer distance to parent coodinate
                          if (dd>weight_limit){
                                        weight_i = 1.0/dd;
                            }else{
                                        weight_i = 1.0/(weight_limit);
                            }
                            weight_tol += weight_i;
                            weight_n++;
                            val_tol += weight_i*T_wall;
                          }
                          
                      }
                  

                      if(weight_n>=1){
                            s[] = val_tol/weight_tol;
                            if(s[]<0){
                                printf("weight_tol=%g\n",weight_tol);
                            }
                      }else{
                          s[] = Tsat0;
                      }
                   }

                }else{
                     s[] =  Tsat0;
                }


          // }//!too_close
       } // !!flag_child_0
       


     }//foreach_child


   }

 
}

attribute {
   bool flag_Tsat0_0;
}
double Tl_refine33_T(Point point, scalar s){  //for T_l

    // write for ff3 need to change when add leafs, so ff3 should be ff only;
    //beside I should write an refine and restriction for centroid_p_l and centroid_p_g
    bool flag_Tsat0_0 = s.flag_Tsat0_0;
     //double Tsat0 = s.Tsat0;
    double T_wall,Tsat0;
    // if(s.boundary[0] != s.boundary_homogeneous[0])
    //     Tsat0 = Tsat00;
    // else
    //     Tsat0 = 0.0;

    // if(s.boundary[0] != s.boundary_homogeneous[0])
    //     T_wall = s.T_wall;
    // else
    //     T_wall = 0.0;

    if(!flag_Tsat0_0)
        Tsat0 = Tsat00;
    else
        Tsat0 = 0.0;

    if(!flag_Tsat0_0)
        T_wall = s.T_wall;
    else
        T_wall = 0.0;




     scalar ff3 = s.ff3;




     double cc;

     double weight_limit =1e-12;

// if(s.inverse){
//      cc = 1.0-ff3[];  
// }else{
//      cc = ff3[];
// }
     cc = coarse(ff3,0,0);

    // double ss = s[];
     //double T_wall = s.T_wall;
     vector centroid_p = s.centroid_p;

     scalar centroidpx = s.centroidpx;
     scalar centroidpy = s.centroidpy;
     scalar centroidpz = s.centroidpz;

     bool flag_wallT = false;
     if(isboundary_right==1 && is_boundary(neighbor(1,0))){
            flag_wallT = true;
     }
     
   if(cc <= 0 ){
     //foreach_child(){
       // s[] = Tsat0; 
        return Tsat0;
    // }
   }else{
     
    
  //     double Delta_coarse = Delta*2.0;

       //coord n = mycs(point,ff3);  // this normal could be changed to being caculate by height later
  //     coord n = interface_normal3(point,ff3,hhh); //normal got from height function.
// if(s.inverse){
//     foreach_dimension(){
//        n.x = -n.x;
//     }
// }
   //    double alpha = plane_alpha (cc, n);

      //for alpha_refine
  //      coord m;
   //     double alphac = 2.*alpha;
   //     foreach_dimension()
    //          m.x = n.x;

 //    foreach_child(){
      //  if((coarse(ff3) > 0. && coarse(ff3,child.x) > 0. &&
      //			   coarse(ff3,0,child.y) > 0. && coarse(ff3,child.x,child.y) > 0.)){
            
      //   }
         // 2-dimensions
        
        coord p_c_child_c; // tansfer to parent's coordate

        // static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
        //     coord nc;
        //     foreach_dimension(){
        //          nc.x = child.x * n.x;
        //     }
        //     double cc2;
        //     cc2 = rectangle_fraction(nc,alpha,a,b);  //child_ff3
           double cc2 = ff3[];
            bool flag_child_0 = false;
          // double alpha_child=alphac;
            double alpha_child;
            coord child_normal;
      if(cc2>EPS && cc2<1.0-EPS){  //if centroid_p refined before s, then there is no need to recaculate p_c_child_c
            //double alpha_child=alphac;
           // foreach_dimension()
                    //alpha_child -= child.x*m.x/2.;  //child_alpha
        
            coord child_p;
            child_normal = interface_normal3(point,ff3,hhh);
            plane_alpha (cc2, child_normal);
            plane_center(child_normal,alpha_child,cc2,&child_p);
           // foreach_dimension(){
                     p_c_child_c.x = child_p.x/2.0 + child.x*0.25;
                     p_c_child_c.y = child_p.y/2.0 + child.y*0.25;
                     p_c_child_c.z = child_p.z/2.0 + child.z*0.25;
           // }
       }else if(cc2>=1){
           p_c_child_c.x = child.x*0.25;
           p_c_child_c.y = child.y*0.25;
           p_c_child_c.z = child.z*0.25;
       }else{ //cc2<=0
         //  s[] = Tsat0;
           flag_child_0 = true;
           return Tsat0;
       }

  if(!flag_child_0){
     double weight_i;
     int ilim = 2*child.x;
     int jlim = 2*child.y;
     int klim = 2*child.z;

      double weight_tol=0.0;
     int weight_n=0;
     double val_tol=0.0;
     bool too_close1=false;
     for(int i=0;abs(i)<abs(ilim);i+=child.x){
        for(int j=0;abs(j)<abs(jlim);j+=child.y){
            for(int k=0;abs(k)<abs(klim);k+=child.z){
              if(!too_close1){         
                        // bool flag_i = true;
                        bool skip1=false;
                        double ssss;
                        //     if(!s.inverse){
                        //           /*
                        // if(coarse(centroid_p.x,i,j)>HUGE/2.0){
                        //     flag_i = false ;
                        //   }
                        // if(coarse(centroid_p.y,i,j)>HUGE/2.0){
                        //     flag_i = false ;
                        // }
                        //           */
                        //         ssss= coarse(ff3,i,j);  
                        //     }else{
                        //         /*
                        //             if(fabs(coarse(centroid_p.x,i,j)-0.0)<1e-20){
                        //     flag_i = false ;
                        //   }
                        // if(fabs(coarse(centroid_p.y,i,j)-0.0)<1e-20){
                        //     flag_i = false ;
                        // }
                        //         */
                        //         ssss = 1.0 - coarse(ff3,i,j); 
                        //     }m
                        ssss= coarse(ff3,i,j,k);
                            if(ssss<=EPS*10){
                                skip1 = true;        
                            }
                            
                      if(!skip1){
                          coord p_temp;
                        //   if(!s.inverse){  
                        //     // double dd = sqrt(sq(p_c_child_c.x-coarse(centroid_p.x,i,j)) + sq(p_c_child_c.y-coarse(centroid_p.y,i,j)));
                        //           //p_temp.x = coarse(centroid_p.x,i,j);
                        //           //p_temp.y = coarse(centroid_p.y,i,j);
                        //           p_temp.x = coarse(centroidpx,i,j);
                        //           p_temp.y = coarse(centroidpy,i,j);
                        //   }else{
                        //     //if(coarse(centroid_p.x,i,j)>HUGE/2.0){
                        //     if(ssss>=1.0){
                        //           p_temp.x = 0.0;
                        //           p_temp.y = 0.0;
                        //     }else{
                        // double sss = coarse(ff3,i,j);
                        //             if(1.0-sss < EPS){
                        //                     // skip1 = true;
                        //                     p_temp.x = 0.0;
                        //                     p_temp.y = 0.0;
                        //             }else{
                        //                     //p_temp.x = -sss*coarse(centroid_p.x,i,j)/(1.0-sss);
                        //   //p_temp.y = -sss*coarse(centroid_p.y,i,j)/(1.0-sss);
                        //                     p_temp.x = -sss*coarse(centroidpx,i,j)/(1.0-sss);
                        //   p_temp.y = -sss*coarse(centroidpy,i,j)/(1.0-sss);
                        //                     foreach_dimension(){
                        //                           if(p_temp.x>=0.5-EPS || p_temp.x<=-0.5+EPS){
                        //                             skip1 = true;
                        //                           }
                        //                     }
                        //             }
                       
                        //     }
                        //   }
                      // if(!s.inverse){  
                                   p_temp.x = coarse(centroidpx,i,j,k);
                                   p_temp.y = coarse(centroidpy,i,j,k);
                                   p_temp.z = coarse(centroidpz,i,j,k);
                      //  }
                        
                            if(!skip1){
                                p_temp.x = i + p_temp.x;
                                p_temp.y = j + p_temp.y;
                                p_temp.z = k + p_temp.z;
                                //double dd = sqrt(sq(p_c_child_c.x-coarse(centroid_p.x,i,j)) + sq(p_c_child_c.y-coarse(centroid_p.y,i,j)));
                                double dd = sqrt(sq(p_c_child_c.x-p_temp.x) + sq(p_c_child_c.y-p_temp.y) + sq(p_c_child_c.z-p_temp.z));
                                if (dd>weight_limit){
                                        weight_i = 1.0/dd;
                                     // if(ssss>EPS){
                                            weight_tol += weight_i;
                                            weight_n++;
                                            val_tol += weight_i*coarse(s,i,j,k);  //20221003
                                     // }
                                }else{
                                        
                                      //  if(ssss>EPS){
                                                too_close1 = true;
                                             //   s[] = coarse(s,i,j);
                                                double temp10=coarse(s,i,j,k);
                                                return temp10;
                                      //  }
                                        //else{
                                         //    s[] = Tsat0*cc2;
                                        //}
                                        //weight_i = 1.0/(weight_limit);
                                }
                               
                            }// compelete1
                      }
              }//too_close
            } //z
        } //y
      }//for child. //x

      // interface_affect
 




 bool interface_effect = true;


         // if(!too_close1){
                bool too_close_to_interface = false;
                //if(cc2>0.0 && cc2<1.0){

 if(interface_effect){
                //if(cc2>=EPS && (cc2<=1-EPS) ){
                if(cc2>=CCCCMIN && (cc2<=1-CCCCMIN) ){
                    coord temp_p = (coord){0.0,0.0,0.0};
                    // double dd = fabs(point_to_line(n,alpha_child,p_c_child_c)); //non-dimensional
                     double dd = fabs(point_to_line(child_normal,alpha_child,temp_p)); //non-dimensional
                    //double dd = fabs(point_to_line(n,alpha,p_c_child_c)); //non-dimensional
                   dd = dd/2.0; // transfer distance to parent coodinate
                    if (dd>weight_limit){
                                weight_i = 1.0/dd;
                    }else{
                                //weight_i = 1.0/(weight_limit);
                                too_close_to_interface = true;
                                
                    }
                    if(!too_close_to_interface){
                        weight_tol += weight_i;
                        weight_n++;
                      // val_tol += weight_i*coarse(s);
                        val_tol += weight_i*Tsat0;
                    }
                }
  }        
                if((!too_close_to_interface)){
                   if(!too_close1){
                      if(flag_wallT){ // right_temperature_wall
                          if(child.y==-1){
                          double dd = fabs(0.5-coarse(centroid_p.x));
                          // dd = dd/2.0; // transfer distance to parent coodinate
                          if (dd>weight_limit){
                                        weight_i = 1.0/dd;
                            }else{
                                        weight_i = 1.0/(weight_limit);
                            }
                            weight_tol += weight_i;
                            weight_n++;
                            val_tol += weight_i*T_wall;
                          }
                          
                      }
                  

                      if(weight_n>=1){
                        // s[] = val_tol/weight_tol*ff3[];
                          // if(weight_tol>weight_limit)
                                //  s[] = val_tol/weight_tol;
                                  return val_tol/weight_tol;
                          // else
                          //       s[] = Tsat0*cc2;

                            if(s[]<0)
                                printf("weight_tol=%g\n",weight_tol);
                      }else{
                        //  s[] = Tsat0*ff3[];
                        //  s[] = Tsat0;
                          return Tsat0;
                        // s[] = Tsat0*0;
                      }
                   }

                }else{
                  //   s[] =  Tsat0;
                     return Tsat0;
                }
          // }//!too_close
       } // !!flag_child_0

          
 //    }//foreach_child


   }

 
}

void Tl_refine3(Point point, scalar s){  //for T_l

    // write for ff3 need to change when add leafs, so ff3 should be ff only;
    //beside I should write an refine and restriction for centroid_p_l and centroid_p_g
     double Tsat0 = s.Tsat0;
     scalar ff3 = s.ff3;

     double cc;

     double weight_limit =1e-12;

if(s.inverse){
     cc = 1.0-ff3[];  
}else{
     cc = ff3[];
}

     double ss = s[];
     double T_wall = s.T_wall;
     vector centroid_p = s.centroid_p;

     scalar centroidpx = s.centroidpx;
     scalar centroidpy = s.centroidpy;

     bool flag_wallT = false;
     if(isboundary_right==1 && is_boundary(neighbor(1,0))){
            flag_wallT = true;
     }
     
   if(cc <= 0 ){
     foreach_child(){
        s[] = Tsat0*0.0; 
     }
   }else{
     
    
     double Delta_coarse = Delta;

       //coord n = mycs(point,ff3);  // this normal could be changed to being caculate by height later
       coord n = interface_normal3(point,ff3,hhh); //normal got from height function.
if(s.inverse){
    foreach_dimension(){
       n.x = -n.x;
    }
}
       double alpha = plane_alpha (cc, n);

      //for alpha_refine
        coord m;
        double alphac = 2.*alpha;
        foreach_dimension()
              m.x = n.x;

     foreach_child(){
      //  if((coarse(ff3) > 0. && coarse(ff3,child.x) > 0. &&
      //			   coarse(ff3,0,child.y) > 0. && coarse(ff3,child.x,child.y) > 0.)){
            
      //   }
         // 2-dimensions
        
        coord p_c_child_c; // tansfer to parent's coordate

        static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
            coord nc;
            foreach_dimension(){
                 nc.x = child.x * n.x;
            }
            double cc2;
            cc2 = rectangle_fraction(nc,alpha,a,b);  //child_ff3
            bool flag_child_0 = false;
           double alpha_child=alphac;
      if(cc2>0.0 && cc2<1.0){  //if centroid_p refined before s, then there is no need to recaculate p_c_child_c
            //double alpha_child=alphac;
            foreach_dimension()
                    alpha_child -= child.x*m.x/2.;  //child_alpha
        
            coord child_p;
            coord child_normal = n;
            plane_center(child_normal,alpha_child,cc2,&child_p);
           // foreach_dimension(){
                     p_c_child_c.x = child_p.x/2.0 + child.x*0.25;
                     p_c_child_c.y = child_p.y/2.0 + child.y*0.25;
           // }
       }else if(cc2>=1){
           p_c_child_c.x = child.x*0.25;
           p_c_child_c.y = child.y*0.25;
       }else{ //cc2<=0
           s[] = Tsat0*0.0;
           flag_child_0 = true;
       }

  if(!flag_child_0){
     double weight_i;
     int ilim = 2*child.x;
     int jlim = 2*child.y;

      double weight_tol=0.0;
     int weight_n=0;
     double val_tol=0.0;
     for(int i=0;abs(i)<abs(ilim);i+=child.x){
        for(int j=0;abs(j)<abs(jlim);j+=child.y){
          // bool flag_i = true;
           bool skip1=false;
           double ssss;
              if(!s.inverse){
                     /*
		       if(coarse(centroid_p.x,i,j)>HUGE/2.0){
		           flag_i = false ;
		        }
		       if(coarse(centroid_p.y,i,j)>HUGE/2.0){
		           flag_i = false ;
		       }
                     */
                  ssss= coarse(ff3,i,j);  
              }else{
                   /*
                      if(fabs(coarse(centroid_p.x,i,j)-0.0)<1e-20){
		           flag_i = false ;
		        }
		      if(fabs(coarse(centroid_p.y,i,j)-0.0)<1e-20){
		           flag_i = false ;
		       }
                   */
                  ssss = 1.0 - coarse(ff3,i,j); 
              }
              if(ssss<=EPS){
                   skip1 = true;        
               }
               
         if(!skip1){
            coord p_temp;
            if(!s.inverse){  
               // double dd = sqrt(sq(p_c_child_c.x-coarse(centroid_p.x,i,j)) + sq(p_c_child_c.y-coarse(centroid_p.y,i,j)));
                    //p_temp.x = coarse(centroid_p.x,i,j);
                    //p_temp.y = coarse(centroid_p.y,i,j);
                    p_temp.x = coarse(centroidpx,i,j);
                    p_temp.y = coarse(centroidpy,i,j);
            }else{
               //if(coarse(centroid_p.x,i,j)>HUGE/2.0){
               if(ssss>=1.0){
                    p_temp.x = 0.0;
                    p_temp.y = 0.0;
               }else{
		       double sss = coarse(ff3,i,j);
                       if(1.0-sss < EPS){
                              // skip1 = true;
                              p_temp.x = 0.0;
                              p_temp.y = 0.0;
                       }else{
                               //p_temp.x = -sss*coarse(centroid_p.x,i,j)/(1.0-sss);
			       //p_temp.y = -sss*coarse(centroid_p.y,i,j)/(1.0-sss);
                               p_temp.x = -sss*coarse(centroidpx,i,j)/(1.0-sss);
			       p_temp.y = -sss*coarse(centroidpy,i,j)/(1.0-sss);
                               foreach_dimension(){
                                    if(p_temp.x>=0.5-EPS || p_temp.x<=-0.5+EPS){
                                       skip1 = true;
                                    }
                               }
                       }
                       
               }
            }

           if(!skip1){
               p_temp.x = i + p_temp.x;
               p_temp.y = j + p_temp.y;
               //double dd = sqrt(sq(p_c_child_c.x-coarse(centroid_p.x,i,j)) + sq(p_c_child_c.y-coarse(centroid_p.y,i,j)));
               double dd = sqrt(sq(p_c_child_c.x-p_temp.x) + sq(p_c_child_c.y-p_temp.y));
               if (dd>weight_limit){
                      weight_i = 1.0/dd;
               }else{
                      weight_i = 1.0/(weight_limit);
               }
               if(ssss>EPS){
                    weight_tol += weight_i;
                    weight_n++;
                    val_tol += weight_i*coarse(s,i,j)/ssss;  //20221003
               }
           }// compelete1
         }
        }
      }//for child.

      // interface_affect
 








    
      //if(cc2>0.0 && cc2<1.0){
      if(cc2>=EPS && (cc2<=1-EPS)){
           double dd = fabs(point_to_line(n,alpha_child,p_c_child_c)); //non-dimensional
          //double dd = fabs(point_to_line(n,alpha,p_c_child_c)); //non-dimensional
          dd = dd/2.0; // transfer distance to parent coodinate
          if (dd>weight_limit){
                      weight_i = 1.0/dd;
          }else{
                      weight_i = 1.0/(weight_limit);
          }
          weight_tol += weight_i;
          weight_n++;
         // val_tol += weight_i*coarse(s);
          val_tol += weight_i*Tsat0;
      }
     
     
     if(flag_wallT){ // right_temperature_wall
        if(child.y==-1){
         double dd = fabs(0.5-coarse(centroid_p.x));
        // dd = dd/2.0; // transfer distance to parent coodinate
         if (dd>weight_limit){
                      weight_i = 1.0/dd;
          }else{
                      weight_i = 1.0/(weight_limit);
          }
          weight_tol += weight_i;
          weight_n++;
          val_tol += weight_i*T_wall;
         }
        
     }
     

         if(weight_n>=1){
           // s[] = val_tol/weight_tol*ff3[];
             // if(weight_tol>weight_limit)
                    s[] = val_tol/weight_tol*cc2;
             // else
             //       s[] = Tsat0*cc2;

              if(s[]<0)
                  printf("weight_tol=%g\n",weight_tol);
         }else{
           //  s[] = Tsat0*ff3[];
            s[] = Tsat0*cc2;
           // s[] = Tsat0*0;
         }
       } // !!flag_child_0

          
     }//foreach_child


   }

 
}

void Tl_restriction3(Point point, scalar s){
//restriction no interface effecct
  scalar ff3 = s.ff3;
  //double cc = ff3[];
  double sum = 0.;
  double cc4;
  double Tsat0 = s.Tsat0;
  if(!s.inverse){
    cc4 = ff3[];
  }else{
    cc4 = 1.0-ff3[];
  }
 // int part_cell=0;
    double weight_i=0.0;
    foreach_child(){
            double cc3;
            if(!s.inverse){
                cc3 = ff3[]; 
            }else{
                cc3 = 1.0 -ff3[];
            }

	    //if(cc3>EPS){
	     //part_cell++;
	    // sum += cc3*s[];  //since
             //sum += s[];
             sum += (s[]*cc3);
             weight_i += cc3 ; 
	    //}
    }
     if(weight_i<=EPS){
         s[] = Tsat0;
     }else{
	  //s[] = sum/weight_i *cc4;
  //  s[] = sum/(1<<dimension);
        s[] = sum/weight_i;
     //s[] = sum*100.0;
     }

     if(s[]<0){

      printf("restriction<0 s[]=%g\n",s[]);
     }
}

void Tl_restriction3_adap(Point point, scalar s){
//restriction no interface effecct
  scalar ff3 = s.ff3;
  //double cc = ff3[];
  double sum = 0.;
  double cc4;
  double Tsat0 = s.Tsat0;
  if(!s.inverse){
    cc4 = ff3[];
  }else{
    cc4 = 1.0-ff3[];
  }
 // int part_cell=0;
    double weight_i=0.0;
    foreach_child(){
            double cc3;
            if(!s.inverse){
                cc3 = ff3[]; 
            }else{
                cc3 = 1.0 -ff3[];
            }

	    //if(cc3>EPS){
	     //part_cell++;
	    // sum += cc3*s[];  //since
             sum += s[]*cc3;
             weight_i += cc3 ; 
	    //}
    }
     if(weight_i<=EPS){
         s[] = Tsat0 ;
     }else{
    //Tl_restriction
   // s[] = sum/(1<<dimension);
      s[] = sum/weight_i;
     }

     if(s[]<0){

      printf("restriction<0 s[]=%g\n",s[]);
     }
}



void centroid_p_x_refine(Point point,scalar s){
     scalar ff3 = s.ff3;
     double cc = ff3[];

if(cc<=0.0 || cc>=1.0){
      foreach_child(){
          s[] = 0.0;
      }
}else{
     coord n = interface_normal3 (point, ff3,hhh);
     
     double alpha = plane_alpha (cc, n);
     coord m;
     double alphac = 2.*alpha;
     foreach_dimension()
              m.x = n.x;
     foreach_child(){
      //  if((coarse(ff3) > 0. && coarse(ff3,child.x) > 0. &&
      //			   coarse(ff3,0,child.y) > 0. && coarse(ff3,child.x,child.y) > 0.)){
            
      //   }
         // 2-dimensions
        
        //coord p_c_child_c; // tansfer to parent's coordate
        coord child_p;
        static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
            coord nc;
            foreach_dimension(){
                 nc.x = child.x * n.x;
            }
            double cc2;
            cc2 = rectangle_fraction(nc,alpha,a,b);  //child_ff3
           double alpha_child=alphac;
      if(cc2>0.0 && cc2<1.0){  //if centroid_p refined before s, then there is no need to recaculate p_c_child_c
            //double alpha_child=alphac;
            foreach_dimension()
                    alpha_child -= child.x*m.x/2.;  //child_alpha
        
           // coord child_p;
            coord child_normal = n;
            plane_center(child_normal,alpha_child,cc2,&child_p);
       }else if(cc2>=1){
           child_p.x = 0.0;
           child_p.y = 0.0;
           child_p.z = 0.0;
       }else{ //cc2<=0
           child_p.x = 0.0;
           child_p.y = 0.0;
           child_p.z = 0.0;
       }
      // foreach_dimension(x){
           s[] = child_p.x;
      // }

    }
  }
}


void centroid_p_y_refine(Point point,scalar s){
     scalar ff3 = s.ff3;
     double cc = ff3[];

if(cc<=0.0 || cc>=1.0){
      foreach_child(){
          s[] = 0.0;
      }
}else{
     coord n = interface_normal3 (point, ff3,hhh);
     
     double alpha = plane_alpha (cc, n);
     coord m;
     double alphac = 2.*alpha;
     foreach_dimension()
              m.x = n.x;
     foreach_child(){
      //  if((coarse(ff3) > 0. && coarse(ff3,child.x) > 0. &&
      //			   coarse(ff3,0,child.y) > 0. && coarse(ff3,child.x,child.y) > 0.)){
            
      //   }
         // 2-dimensions
        
        //coord p_c_child_c; // tansfer to parent's coordate
        coord child_p;
        static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
            coord nc;
            foreach_dimension(){
                 nc.x = child.x * n.x;
            }
            double cc2;
            cc2 = rectangle_fraction(nc,alpha,a,b);  //child_ff3
           double alpha_child=alphac;
      if(cc2>0.0 && cc2<1.0){  //if centroid_p refined before s, then there is no need to recaculate p_c_child_c
            //double alpha_child=alphac;
            foreach_dimension()
                    alpha_child -= child.x*m.x/2.;  //child_alpha
        
           // coord child_p;
            coord child_normal = n;
            plane_center(child_normal,alpha_child,cc2,&child_p);
       }else if(cc2>=1){
           child_p.x = 0.0;
           child_p.y = 0.0;
           child_p.z = 0.0;
       }else{ //cc2<=0
           child_p.x = 0.0;
           child_p.y = 0.0;
           child_p.z = 0.0;
       }
     //  foreach_dimension(y){
           s[] = child_p.y;
     //  }

    }
  }
}

void centroid_p_z_refine(Point point,scalar s){
     scalar ff3 = s.ff3;
     double cc = ff3[];

if(cc<=0.0 || cc>=1.0){
      foreach_child(){
          s[] = 0.0;
      }
}else{
     coord n = interface_normal3 (point, ff3,hhh);
     
     double alpha = plane_alpha (cc, n);
     coord m;
     double alphac = 2.*alpha;
     foreach_dimension()
              m.x = n.x;
     foreach_child(){
      //  if((coarse(ff3) > 0. && coarse(ff3,child.x) > 0. &&
      //			   coarse(ff3,0,child.y) > 0. && coarse(ff3,child.x,child.y) > 0.)){
            
      //   }
         // 2-dimensions
        
        //coord p_c_child_c; // tansfer to parent's coordate
        coord child_p;
        static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
            coord nc;
            foreach_dimension(){
                 nc.x = child.x * n.x;
            }
            double cc2;
            cc2 = rectangle_fraction(nc,alpha,a,b);  //child_ff3
           double alpha_child=alphac;
      if(cc2>0.0 && cc2<1.0){  //if centroid_p refined before s, then there is no need to recaculate p_c_child_c
            //double alpha_child=alphac;
            foreach_dimension()
                    alpha_child -= child.x*m.x/2.;  //child_alpha
        
           // coord child_p;
            coord child_normal = n;
            plane_center(child_normal,alpha_child,cc2,&child_p);
       }else if(cc2>=1){
           child_p.x = 0.0;
           child_p.y = 0.0;
           child_p.z = 0.0;
       }else{ //cc2<=0
           child_p.x = 0.0;
           child_p.y = 0.0;
           child_p.z = 0.0;
       }
     //  foreach_dimension(y){
           s[] = child_p.z;
     //  }

    }
  }
}

void centroid_p_x_restriction(Point point,scalar s){
     double val_total = 0.0;
     double v_total=0.0;
     int val_n=0;
     double weight_i;
     double weight_tol=0.0;
     scalar ff3 = s.ff3;
     bool flag4 = true;
     foreach_child(){
        if(ff3[]<1.0 && ff3[]>0){
            flag4 = false;
        }
        if(ff3[]>0.0){
            //if(s[]<0.5 && s[]>-0.5){
                //p_c_child_c.x = child_p.x/2.0 + child.x*0.25;
                //p_c_child_c.y = child_p.y/2.0 + child.y*0.25;
                double temp3 = s[]/2.0+child.x*0.25;
                val_n++;
                weight_i = ff3[];
                weight_tol += weight_i;
                val_total += temp3*weight_i;
                
            // }
        }
     }
     if(flag4){
        s[] = 0.0;
     }else{
        if(val_n>=1){
          // s[] = val_total/weight_tol;
           s[] = val_total/(1<<dimension);
        }else{
           s[] = 0.0;
        }
     }
}

void centroid_p_y_restriction(Point point,scalar s){
     double val_total = 0.0;
     double v_total=0.0;
     int val_n=0;
     double weight_i;
     double weight_tol=0.0;
     scalar ff3 = s.ff3;
     bool flag4 = true;
     foreach_child(){
        if(ff3[]<=1.0 && ff3[]>=0){
            flag4 = false;
        }
        if(ff3[]>0.0){
            //if(s[]<0.5 && s[]>-0.5){
                //p_c_child_c.x = child_p.x/2.0 + child.x*0.25;
                //p_c_child_c.y = child_p.y/2.0 + child.y*0.25;
                double temp3 = s[]/2.0+child.y*0.25;
                val_n++;
                weight_i = ff3[];
                weight_tol += weight_i;
                val_total += temp3*weight_i;
                
            // }
        }
     }
     if(flag4){
        s[] = 0.0;
     }else{
        if(val_n>=1){
           // s[] = val_total/weight_tol;
              s[] = val_total/(1 << dimension);
        }else{
           s[] = 0.0;
        }
     }
}

void centroid_p_z_restriction(Point point,scalar s){
     double val_total = 0.0;
     double v_total=0.0;
     int val_n=0;
     double weight_i;
     double weight_tol=0.0;
     scalar ff3 = s.ff3;
     bool flag4 = true;
     foreach_child(){
        if(ff3[]<=1.0 && ff3[]>=0){
            flag4 = false;
        }
        if(ff3[]>0.0){
            //if(s[]<0.5 && s[]>-0.5){
                //p_c_child_c.x = child_p.x/2.0 + child.x*0.25;
                //p_c_child_c.y = child_p.y/2.0 + child.y*0.25;
                double temp3 = s[]/2.0+child.z*0.25;
                val_n++;
                weight_i = ff3[];
                weight_tol += weight_i;
                val_total += temp3*weight_i;
                
            // }
        }
     }
     if(flag4){
        s[] = 0.0;
     }else{
        if(val_n>=1){
           // s[] = val_total/weight_tol;
              s[] = val_total/(1 << dimension);
        }else{
           s[] = 0.0;
        }
     }
}


extern scalar centroidpx,centroidpy,centroidpz;
extern vector hhh;
void T_tree_function(scalar TT, scalar ff, double Tsat0){


/*
     scalar ff3[];
     if(TT.inverse){
          foreach(){
             ff3[] = 1.0-ff[];
          }
      }else{
          foreach(){
             ff3[] = ff[];
          }
      }
*/

  /*
     if(flag1){
          foreach(){
             ff3[] = 1.0-ff[];
          }
      }else{
          foreach(){
             ff3[] = ff[];
          }
      }
  */
      
     //for (scalar c in interfaces)
       //for (scalar c in {ff3})
	    //if (ff3.height.x.i)
	      //heights (ff3,ff3.height);
        // if (ff.height.x.i)
	 //     heights (ff,ff.height);

	 // for centroid_p
        
         heights (ff,hhh);
	vector centroid_p[];

	  foreach_level(0){
              centroidpx[] = nodata;
              centroidpy[] = nodata;
              centroidpz[] = nodata;
	      foreach_dimension()
		centroid_p.x[] = nodata;
          }
	  
	  //restriction ({ff3});
    //ff.restriction = restriction_volume_average;

          restriction ({ff});
	    for (int l = 0; l <= depth(); l++) {
	       foreach_level(l){
                  //calculate centroid p
              // double cc = ff3[];
                            double cc = ff[];
              if(cc>0.0 && cc<1.0){
                  coord p_c;//p_c_real;
                  //coord n = mycs(point,ff3);  // this normal could be changed to being caculate by height later
                              coord n = interface_normal3(point,ff,hhh);  
                  double alpha = plane_alpha (cc, n);
                  plane_center(n,alpha,cc,&p_c); //center of liquid
                  foreach_dimension(){
                      centroid_p.x[] = p_c.x;
                  }
                              centroidpx[] = p_c.x;
                              centroidpy[] = p_c.y;
                              centroidpz[] = p_c.z;
                }else if(cc>=1){
                      foreach_dimension(){
                          centroid_p.x[] = 0.0;
                      }
                                  centroidpx[] = 0.0;
                                  centroidpy[] = 0.0;
                                  centroidpz[] = 0.0;
                }else{  //cc <=0
                      foreach_dimension(){
                          centroid_p.x[] = 0.0; //HUGE
                      }
                                  centroidpx[] = 0.0;
                                  centroidpy[] = 0.0;
                                  centroidpz[] = 0.0;
                }
		    
	       }
               if(l>=1){
                   foreach_halo(prolongation,l-1){
                       //assume that resolution boundary are all pure cell
                      foreach_child(){
                        foreach_dimension(){
                          centroid_p.x[] = 0.0;
                        }
                          centroidpx[] = 0.0;
                          centroidpy[] = 0.0;
                          centroidpz[] = 0.0;
                      } 
                   }//foreach_halo
               }//for halo
	    }

        restriction({centroidpx,centroidpy,centroidpz,centroid_p.x});

        centroidpx.ff3=ff;
        centroidpy.ff3=ff;
        centroidpz.ff3=ff;
        centroidpx.refine = centroid_p_x_refine;
        centroidpy.refine = centroid_p_y_refine;
        centroidpz.refine = centroid_p_z_refine; 

        centroidpx.restriction = centroid_p_x_restriction;   //I am here
        centroidpy.restriction = centroid_p_y_restriction;
        centroidpz.restriction = centroid_p_z_restriction;
        centroidpx.dirty=true;
        centroidpy.dirty=true;
        centroidpz.dirty=true;

       // restriction({centroidpx,centroidpy,centroid_p.x});
/*         
	TT.ff3 = ff3;
	TT.centroid_p=centroid_p;
	TT.Tsat0 = Tsat0;
	TT.T_wall = boundary_value_right;

	TT.refine = TT.prolongation = Tl_refine2;
	TT.restriction = Tl_restriction2;
        TT.dirty = true;
*/
        TT.hhh = hhh;
        TT.ff3 = ff;
	TT.centroid_p=centroid_p;
        TT.centroidpx=centroidpx;
        TT.centroidpy=centroidpy;
        TT.centroidpz=centroidpz;

	TT.Tsat0 = Tsat0;
	TT.T_wall = boundary_value_right;


	TT.refine = TT.prolongation = Tl_refine33; //Tl_refine33; // changed for Tlff
	TT.restriction = Tl_restriction3;
        TT.dirty = true;
}




void T_tree_function_adap(scalar TT, scalar ff, double Tsat0){

        
  heights (ff,hhh);
	vector centroid_p[];

	  foreach_level(0){
              centroidpx[] = nodata;
              centroidpy[] = nodata;
              centroidpz[] = nodata;
	      foreach_dimension()
		centroid_p.x[] = nodata;
          }
	  
	  //restriction ({ff3});
    //ff.restriction = restriction_volume_average;

          restriction ({ff});
	    for (int l = 0; l <= depth(); l++) {
	       foreach_level(l){
		    //calculate centroid p
		 // double cc = ff3[];
                  double cc = ff[];
		 if(cc>0.0 && cc<1.0){
		     coord p_c;//p_c_real;
		     //coord n = mycs(point,ff3);  // this normal could be changed to being caculate by height later
                     coord n = interface_normal3(point,ff,hhh);  
		     double alpha = plane_alpha (cc, n);
		     plane_center(n,alpha,cc,&p_c); //center of liquid
		     foreach_dimension(){
		         centroid_p.x[] = p_c.x;
		     }
                     centroidpx[] = p_c.x;
                     centroidpy[] = p_c.y;
                     centroidpz[] = p_c.z;
		  }else if(cc>=1){
		         foreach_dimension(){
		             centroid_p.x[] = 0.0;
		         }
                         centroidpx[] = 0.0;
                         centroidpy[] = 0.0;
                         centroidpz[] = 0.0;
		  }else{  //cc <=0
		         foreach_dimension(){
		             centroid_p.x[] = 0.0; //HUGE
		         }
                         centroidpx[] = 0.0;
                         centroidpy[] = 0.0;
                         centroidpz[] = 0.0;
		  }
		    
	       }
               if(l>=1){
                   foreach_halo(prolongation,l-1){
                       //assume that resolution boundary are all pure cell
                      foreach_child(){
		                     foreach_dimension(){
				                    centroid_p.x[] = 0.0;
				              }
		                 centroidpx[] = 0.0;
		                 centroidpy[] = 0.0;
                     centroidpz[] = 0.0;
                      } 
                   }//foreach_halo
               }//for halo
	    }

        restriction({centroidpx,centroidpy,centroidpz,centroid_p.x});

        centroidpx.ff3=ff;
        centroidpy.ff3=ff;
        centroidpz.ff3=ff;
        centroidpx.refine = centroid_p_x_refine;
        centroidpy.refine = centroid_p_y_refine;
        centroidpz.refine = centroid_p_z_refine; 

        centroidpx.restriction = centroid_p_x_restriction;   //I am here
        centroidpy.restriction = centroid_p_y_restriction;
        centroidpz.restriction = centroid_p_z_restriction;
        centroidpx.dirty=true;
        centroidpy.dirty=true;
        centroidpz.dirty=true;


        TT.hhh = hhh;
        TT.ff3 = ff;
	TT.centroid_p=centroid_p;
        TT.centroidpx=centroidpx;
        TT.centroidpy=centroidpy;
        TT.centroidpz=centroidpz;
	TT.Tsat0 = Tsat0;
	TT.T_wall = boundary_value_right;  

	TT.refine = TT.prolongation = Tl_refine33_adap; //Tl_refine33;  // _adap checked for   3D

	TT.restriction = Tl_restriction3_adap;  //checked for 3D
        TT.dirty = true;
}

//this is used in poisson3.h; there ff is 1.0-ff[].
//so correspondingly we should modify Tl_refine33 to Tl_refine33_T
void T_tree_function_T_poisson(scalar TT, scalar ff, double Tsat0){

         heights (ff,hhh);
	vector centroid_p[];

	  foreach_level(0){
              centroidpx[] = nodata;
              centroidpy[] = nodata;
              centroidpz[] = nodata;
	      foreach_dimension()
		centroid_p.x[] = nodata;
          }
	  
	  //restriction ({ff3});
          restriction ({ff});
	    for (int l = 0; l <= depth(); l++) {
	       foreach_level(l){
		    //calculate centroid p
		 // double cc = ff3[];
                  double cc = ff[];
		 if(cc>0.0 && cc<1.0){
		     coord p_c;//p_c_real;
		     //coord n = mycs(point,ff3);  // this normal could be changed to being caculate by height later
                     coord n = interface_normal3(point,ff,hhh);  
		     double alpha = plane_alpha (cc, n);
		     plane_center(n,alpha,cc,&p_c); //center of liquid
		     foreach_dimension(){
		         centroid_p.x[] = p_c.x;
		     }
                     centroidpx[] = p_c.x;
                     centroidpy[] = p_c.y;
                     centroidpz[] = p_c.z;
		  }else if(cc>=1){
		         foreach_dimension(){
		             centroid_p.x[] = 0.0;
		         }
                         centroidpx[] = 0.0;
                         centroidpy[] = 0.0;
                         centroidpz[] = 0.0;
		  }else{  //cc <=0
		         foreach_dimension(){
		             centroid_p.x[] = 0.0; //HUGE
		         }
                         centroidpx[] = 0.0;
                         centroidpy[] = 0.0;
                         centroidpz[] = 0.0;
		  }
		    
	       }
               if(l>=1){
                   foreach_halo(prolongation,l-1){
                       //assume that resolution boundary are all pure cell
                      foreach_child(){
                            foreach_dimension(){
                                    centroid_p.x[] = 0.0;
                            }
                              centroidpx[] = 0.0;
                              centroidpy[] = 0.0;
                              centroidpz[] = 0.0;
                      } 
                   }//foreach_halo
               }//for halo
	    }

      //  restriction({centroidpx,centroidpy,centroid_p.x});

     //   centroidpx.ff3=ff;
     //   centroidpy.ff3=ff;
     //   centroidpx.refine = centroid_p_x_refine;
     //   centroidpy.refine = centroid_p_y_refine;
     //   centroidpx.restriction = centroid_p_x_restriction;
      //  centroidpy.restriction = centroid_p_y_restriction;
     //   centroidpx.dirty=true;
     //   centroidpy.dirty=true;
/*         
	TT.ff3 = ff3;
	TT.centroid_p=centroid_p;
	TT.Tsat0 = Tsat0;
	TT.T_wall = boundary_value_right;

	TT.refine = TT.prolongation = Tl_refine2;
	TT.restriction = Tl_restriction2;
        TT.dirty = true;
*/
    // restriction({centroidpx,centroidpy,ff});

        TT.hhh = hhh;
        TT.ff3 = ff;
	//TT.centroid_p=centroid_p;
        TT.centroidpx=centroidpx;
        TT.centroidpy=centroidpy;
        TT.centroidpz=centroidpz;
      
	//TT.Tsat0 = Tsat0;
	TT.T_wall = boundary_value_right;

  

	//TT.refine = TT.prolongation = Tl_refine33_T;
	//TT.restriction = Tl_restriction3;
   //     TT.dirty = true;
}


//void T_tree_function2(scalar TT, scalar ff){

//      TT.ff2 = ff;
//      TT.Trefine = T.prolongation = energy_refine3 ;
//      TT.restriction = restriction_cprho_average3 ;
     
//}



// need height, restriction(ff),Tsat0,Twall;
void T_refine3(Point point, scalar s){  //for T

     double Tsat0 = s.Tsat0;

     scalar ff3 = s.ff3;

     double cc;
     cc = ff3[];
//  which side T is used to refine is determined by volume of volume of child cell;

     double ss = s[];
     double T_wall = s.T_wall;


     bool flag_wallT = false;
     if(isboundary_right==1 && is_boundary(neighbor(1,0))){
            flag_wallT = true;
     }
     

     
    
     double Delta_coarse = Delta;

       
       coord n = interface_normal3(point,ff3,hhh); //normal got from height function.


      // double alpha = plane_alpha (cc, n);
       
       double alpha = plane_alpha (cc, n);
       coord m;
       double alphac = 2.*alpha;
       foreach_dimension()
              m.x = n.x;

       
     foreach_child(){

       

       double cc2;
       
         // 2-dimensions
       bool phase_flag=true; // true: liquid; false: gas;
       if(cc<=0){ //parent is pure gas, child is gas
            phase_flag=false;      
            cc2 =0.0;
       }else{
         

         static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
            coord nc;
            foreach_dimension(){
                 nc.x = child.x * n.x;
            }
            
            cc2 = rectangle_fraction(nc,alpha,a,b);  //child_ff3
        
        if(cc2<0.5){
            phase_flag = false;
         }
    
      }

           coord p_c_child_c; // tansfer to parent's coordate

        
            coord child_p;
            child_p=(coord){0.0,0.0};

                     p_c_child_c.x = child_p.x/2.0 + child.x*0.25;
                     p_c_child_c.y = child_p.y/2.0 + child.y*0.25;

     double weight_i;
     int ilim = 2*child.x;
     int jlim = 2*child.y;

     double weight_tol=0.0;
     int weight_n=0;
     double val_tol=0.0;
     for(int i=0;abs(i)<abs(ilim);i+=child.x){
        for(int j=0;abs(j)<abs(jlim);j+=child.y){
         
           bool skip1=false;
           double ssss;

            ssss = coarse(ff3,i,j);
            if ( !((ssss < 0.5 && !phase_flag) || (ssss >= 0.5 && phase_flag))){
                 skip1 = true;
            }  
         if(!skip1){
           coord p_temp;
           foreach_dimension(){
               p_temp.x = 0.0;
           }

               p_temp.x = i + p_temp.x; //coarse cell coordinate relative to children's parent
               p_temp.y = j + p_temp.y;
               //double dd = sqrt(sq(p_c_child_c.x-coarse(centroid_p.x,i,j)) + sq(p_c_child_c.y-coarse(centroid_p.y,i,j)));
               double dd = sqrt(sq(p_c_child_c.x-p_temp.x) + sq(p_c_child_c.y-p_temp.y));
               if (dd>1e-12){
                      weight_i = 1.0/dd;
               }else{
                      weight_i = 1.0/(1e-12);
               }
               weight_tol += weight_i;
               weight_n++;
               val_tol += weight_i*coarse(s,i,j);
         }
        }
      }//for child.

      // interface_affect

    
    if(cc2>0.0 && cc2<1.0){
           double alpha_child=alphac;
           foreach_dimension()
                    alpha_child -= child.x*m.x/2.;  //child_alpha
           double dd = fabs(point_to_line(n,alpha_child,p_c_child_c)); //non-dimensional
          //double dd = fabs(point_to_line(n,alpha,p_c_child_c)); //non-dimensional
          dd = dd/2.0; // transfer distance to parent coodinate
          if (dd>1e-12){
                      weight_i = 1.0/dd;
          }else{
                      weight_i = 1.0/(1e-12);
          }
          weight_tol += weight_i;
          weight_n++;
          val_tol += weight_i*coarse(s);
      }
     
     
     if(flag_wallT){ // right_temperature_wall
        if(child.y==-1){
         double dd = fabs(0.5-0.0);
         dd = dd/2.0; // transfer distance to parent coodinate
         if (dd>1e-12){
                      weight_i = 1.0/dd;
          }else{
                      weight_i = 1.0/(1e-12);
          }
          weight_tol += weight_i;
          weight_n++;
          val_tol += weight_i*T_wall;
         }
        
     }

         if(weight_n>=1){
           // s[] = val_tol/weight_tol*ff3[];
             s[] = val_tol/weight_tol;
         }else{
           //  s[] = Tsat0*ff3[];
           s[] = Tsat0;  //????????????
         }
     //  } // !!flag_child_0

          
     }//foreach_child


  // }

}


/*
//task for 20220914
// modify T_refine3 to T_refine3_0  and T_refine3_1. Difference is that when cc<0 or cc2<0.5, T=Tsat0.  FINISHED
// add T_restriction3(_0,T_resriction3_1). note, when cc<0.5, T=Tsat0;                                  FINISHED
// solve T two time, one for T gas and another is T liquid. so need T_old to save original value        

attribute {
    bool poisson_phase; // true: fluid side; false: gas side 
}
// need height, restriction(ff),poisson_phase,Tsat0,Twall;
void T_refine3_2(Point point, scalar s){  //for T refine in poisson solver

    double Tsat0 = s.Tsat0;
  
     bool poisson_phase = s.poisson_phase;
     scalar ff3 = s.ff3;

     double cc;
     if(!poisson_phase){
	  cc = 1.0 - ff3[];
     }else{
          cc = ff3[];
     }
//  which side T is used to refine is determined by volume of volume of child cell;

     double ss = s[];
     double T_wall = s.T_wall;



     bool flag_wallT = false;
     if(isboundary_right==1 && is_boundary(neighbor(1,0))){
            flag_wallT = true;
     }
     

     
     
     double Delta_coarse = Delta;

       
       coord n = interface_normal3(point,ff3,hhh); //normal got from height function.


       double alpha = plane_alpha (cc, n);



       coord m;
       double alphac = 2.*alpha;
       foreach_dimension()
              m.x = n.x;

     foreach_child(){
       double cc2;
       bool complete1=false;
         // 2-dimensions
       bool phase_flag=poisson_phase; //initial true: liquid; false: gas;
       if(cc<=0){ //parent is pure gas, child is gas
            phase_flag=!poisson_phase;
            s[] = Tsat0;
            complete1 = true;      
       }else{
         

         static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
            coord nc;
            foreach_dimension(){
                 nc.x = child.x * n.x;
            }
            
            cc2 = rectangle_fraction(nc,alpha,a,b);  //child_ff3
            if(!poisson_phase){
                cc2 = 1.0 - cc2;
            } 
            if(cc2<0.5){
              phase_flag = !poisson_phase;
              s[] = Tsat0;
              complete1 = true;
            }
    
       }


      if(!complete1){
            coord p_c_child_c; // tansfer to parent's coordate
            coord child_p;
            child_p=(coord){0.0,0.0};

                     p_c_child_c.x = child_p.x/2.0 + child.x*0.25;
                     p_c_child_c.y = child_p.y/2.0 + child.y*0.25;

            double weight_i;
            int ilim = 2*child.x;
            int jlim = 2*child.y;

            double weight_tol=0.0;
            int weight_n=0;
            double val_tol=0.0;
	     for(int i=0;abs(i)<abs(ilim);i+=child.x){
		for(int j=0;abs(j)<abs(jlim);j+=child.y){
		 
		   bool skip1=false;
		   double ssss;

		    ssss = coarse(ff3,i,j);
                    if(!phase_flag){
                        ssss = 1.0 - ssss;
                    } 
                    if(ssss<0.5){
                         skip1 = true;
                     } 
		 if(!skip1){
		     coord p_temp;
		     foreach_dimension(){
		         p_temp.x = 0.0;
		     }

		       p_temp.x = i + p_temp.x;
		       p_temp.y = j + p_temp.y;
		       //double dd = sqrt(sq(p_c_child_c.x-coarse(centroid_p.x,i,j)) + sq(p_c_child_c.y-coarse(centroid_p.y,i,j)));
		       double dd = sqrt(sq(p_c_child_c.x-p_temp.x) + sq(p_c_child_c.y-p_temp.y));
		       if (dd>1e-12){
		              weight_i = 1.0/dd;
		       }else{
		              weight_i = 1.0/(1e-12);
		       }
		       weight_tol += weight_i;
		       weight_n++;
		       val_tol += weight_i*coarse(s,i,j);
		 }
		}
	      }//for child.

	      // interface_affect
	    if(cc2>0.0 && cc2<1.0){  // actrually cc2 must >0.5
                  double alpha_child=alphac;
                  foreach_dimension()
                    alpha_child -= child.x*m.x/2.;  //child_alpha
		  // double dd = fabs(point_to_line(child_normal,alpha_child,p_c_child_c)); //non-dimensional
		  double dd = fabs(point_to_line(n,alpha,p_c_child_c)); //non-dimensional
		  //dd = dd/2.0; // transfer distance to parent coodinate
		  if (dd>1e-12){
		              weight_i = 1.0/dd;
		  }else{
		              weight_i = 1.0/(1e-12);
		  }
		  weight_tol += weight_i;
		  weight_n++;
		  val_tol += weight_i*coarse(s);
	      }
	     
	     
	     if(flag_wallT){ // right_temperature_wall
		if(child.y==-1){
		 double dd = fabs(0.5-0.0);
		 dd = dd/2.0; // transfer distance to parent coodinate
		 if (dd>1e-12){
		              weight_i = 1.0/dd;
		  }else{
		              weight_i = 1.0/(1e-12);
		  }
		  weight_tol += weight_i;
		  weight_n++;
		  val_tol += weight_i*T_wall;
		 }
		
	     }

		 if(weight_n>=1){
		   // s[] = val_tol/weight_tol*ff3[];
		     s[] = val_tol/weight_tol;
		 }else{
		   //  s[] = Tsat0*ff3[];
		   s[] = Tsat0;  //????????????
		 }
	     //  } // !!flag_child_0

         }
     }//foreach_child


  // }

}


//restriction function in poisson solver
//before using this, remember to restriction ff
// need height, restriction(ff),poisson_phase,Tsat0;
void T_restriction3_2(Point point, scalar s){
  bool poisson_phase = s.poisson_phase;
  scalar ff3 = s.ff3;
  double Tsat0 = s.Tsat0;
  //double cc = ff3[];
  double sum = 0.;
  double cc4;
  bool complete1=false;
  if(poisson_phase){
    cc4 = ff3[];
  }else{
    cc4 = 1.0-ff3[];
  }
  if(cc4<0.5){
      s[] = Tsat0;
      complete1 = true;
   }
  
  if(!complete1){
        
        
        
        int weight_n = 0;
	foreach_child(){
             double cc3;
             if(poisson_phase){
		        cc3 = ff3[]; 
	      }else{
		        cc3 = 1.0 -ff3[];
	      }
	      if(cc3>0){
		     //part_cell++;
		 // sum += cc3*s[];  //since
                 double dd = sqrt(sq(child.x*0.25)+sq(child.y*0.25));
                
                 sum += 1.0/dd*s[];
                 weight_n++; 
	       }
	 }
         if(weight_n<pow(2.0,dimension)){ //if <4 for 2-dimension/ calculate interface efact
         //only need to contain the effect of interface; because Twall interface is on one side, there will 
         // always be children or interface
            coord n = interface_normal3 (point, ff3,hhh),p;
            double alpha = plane_alpha (ff3[], n);
            plane_area_center(n,alpha,&p);
            double dd = sqrt(sq(p.x)+sq(p.y));
            if (dd<1e-13){
                   dd = 1e-13;
            }
            sum += s[]*1.0/dd;
            weight_n++;
         }
         if(weight_n>=1){
              s[] = sum/(1 << dimension)/(cc4 + 1e-30)*cc4;
              complete1 = true;
         }else{
                printf(stderr,"stop restriction error!!!!!!!\n");
         }
   }

}

//use this in poisson solver

void T_tree_function_for_T(scalar TT, scalar ff, bool poisson_phase){

        // if (ff.height.x.i)
	//      heights (ff,ff.height);
       
        
        
        heights (ff,hhh);
        restriction({ff});
        TT.hhh=hhh;
        TT.poisson_phase = poisson_phase;
        TT.ff3 = ff;
	TT.Tsat0 = Tsat0;
	TT.T_wall = boundary_value_right;

	TT.refine = TT.prolongation = T_refine3_2;
	TT.restriction = T_restriction3_2;
        TT.dirty = true;
}
*/
