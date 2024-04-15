/**
# Contact angles on an embedded boundary

This header file implemented contact angles for [VOF
interfaces](/src/vof.h) on [embedded boundaries](/src/embed.h).

The contact angle is defined by this field and can be constant or
variable. */

#include "./reconstruction4.h" //reconstruct using height_normal

extern scalar intersect_true;
extern int level_interface;

extern scalar corner_ff;

(const) scalar contact_angle;
//scalar contact_angle[];

/**
This function returns the properly-oriented normal of an interface
touching the embedded boundary. `ns` is the normal to the embedded
boundary, `nf` the VOF normal, not taking into account the contact
angle, and `angle` the contact angle. */

// static inline coord normal_contact (coord ns, coord nf, double angle)
// {
//   coord n;
//   if (dimension == 2){
//     if (- ns.x*nf.y + ns.y*nf.x > 0) { // 2D cross product
//       n.x = - ns.x*cos(angle) + ns.y*sin(angle);
//       n.y = - ns.x*sin(angle) - ns.y*cos(angle);
//     }
//     else {
//       n.x = - ns.x*cos(angle) - ns.y*sin(angle);
//       n.y =   ns.x*sin(angle) - ns.y*cos(angle);
//     }
//   }
//   else{
//     // fixme: Only horizontal embeld plane for the moment 
//     //if(fabs(nf.z)<fabs(nf.x)){
//     if(1==1){
//         double phi = atan2(nf.x,nf.z);
//         n.x = sin(angle)*sin(phi) ;
//         n.y = cos(angle) ;
//         n.z = sin(angle)*cos(phi);
//     }else{
//         double phi = atan2(nf.z,nf.x);
//         n.x = sin(angle)*cos(phi) ;
//         n.y = cos(angle) ;
//         n.z = sin(angle)*sin(phi);
//     }

//      n = nf;
//   }
//   return n;

// }
// // static inline coord normal_contact (coord ns, coord nf, double angle)
// // {
// //   coord nc; 
// //   #if dimension == 2
// //   if (- ns.x*nf.y + ns.y*nf.x > 0) { // 2D cross product // difference between pointing towards or away
    
// //     //  standard matrix *-1
// //     nc.x = - ns.x*cos(angle) + ns.y*sin(angle);
// //     nc.y = - ns.x*sin(angle) - ns.y*cos(angle);
// //   }
// //   else { // inverted
// //     nc.x = - ns.x*cos(angle) - ns.y*sin(angle);
// //     nc.y = + ns.x*sin(angle) - ns.y*cos(angle);
// //   }
 
// //   #elif dimension == 3
// //   if(1==0){
// //           coord rot_axis;
// //           // normal vector of ns and nf calculated by cross product.
// //           rot_axis.x = -ns.y*nf.z + ns.z*nf.y;
// //           rot_axis.y = -ns.z*nf.x + ns.x*nf.z;
// //           rot_axis.z = -ns.x*nf.y + ns.y*nf.x;
// //           normalize(&rot_axis);
// //           // rotate ns by angle around rot_axis
// //           nc.x = - ns.x*(sq(rot_axis.x)          *(1.-cos(angle)) +             cos(angle)) 
// //                 - ns.y*(rot_axis.x*rot_axis.y*(1.-cos(angle)) - rot_axis.z*sin(angle)) 
// //                 - ns.z*(rot_axis.x*rot_axis.z*(1.-cos(angle)) + rot_axis.y*sin(angle));
              
// //           nc.y = - ns.x*(rot_axis.y*rot_axis.x*(1.-cos(angle)) + rot_axis.z*sin(angle)) 
// //                 - ns.y*(sq(rot_axis.y)       *(1.-cos(angle)) +            cos(angle)) 
// //                 - ns.z*(rot_axis.y*rot_axis.z*(1.-cos(angle)) -rot_axis.x*sin(angle));
              
// //           nc.z = - ns.x*(rot_axis.z*rot_axis.x*(1.-cos(angle)) - rot_axis.y*sin(angle)) 
// //                 - ns.y*(rot_axis.z*rot_axis.y*(1.-cos(angle)) + rot_axis.x*sin(angle)) 
// //                 - ns.z*(sq(rot_axis.z)       *(1.-cos(angle)) +            cos(angle));  
// //   }else{
// //           coord rot_axis;
// //           coord mid_vector;
// //           double val1=0.0;;
// //           foreach_dimension(){
// //             val1=val1+nf.x*ns.x;
// //           }
// //           //mid_vecor nf - (nf*ns)*ns
// //           foreach_dimension(){
// //               mid_vector.x = nf.x - val1*ns.x;
// //           }
// //           normalize(&mid_vector);
// //           // normal vector of ns and nf calculated by cross product.
// //           // rot_axis.x = -ns.y*nf.z + ns.z*nf.y;
// //           // rot_axis.y = -ns.z*nf.x + ns.x*nf.z;
// //           // rot_axis.z = -ns.x*nf.y + ns.y*nf.x;
// //           rot_axis.x = -ns.y*mid_vector.z + ns.z*mid_vector.y;
// //           rot_axis.y = -ns.z*mid_vector.x + ns.x*mid_vector.z;
// //           rot_axis.z = -ns.x*mid_vector.y + ns.y*mid_vector.x;
// //           normalize(&rot_axis);
// //           // rotate ns by angle around rot_axis
// //           nc.x = - ns.x*(sq(rot_axis.x)          *(1.-cos(angle)) +             cos(angle)) 
// //                 - ns.y*(rot_axis.x*rot_axis.y*(1.-cos(angle)) - rot_axis.z*sin(angle)) 
// //                 - ns.z*(rot_axis.x*rot_axis.z*(1.-cos(angle)) + rot_axis.y*sin(angle));
              
// //           nc.y = - ns.x*(rot_axis.y*rot_axis.x*(1.-cos(angle)) + rot_axis.z*sin(angle)) 
// //                 - ns.y*(sq(rot_axis.y)       *(1.-cos(angle)) +            cos(angle)) 
// //                 - ns.z*(rot_axis.y*rot_axis.z*(1.-cos(angle)) -rot_axis.x*sin(angle));
              
// //           nc.z = - ns.x*(rot_axis.z*rot_axis.x*(1.-cos(angle)) - rot_axis.y*sin(angle)) 
// //                 - ns.y*(rot_axis.z*rot_axis.y*(1.-cos(angle)) + rot_axis.x*sin(angle)) 
// //                 - ns.z*(sq(rot_axis.z)       *(1.-cos(angle)) +            cos(angle));  
// //   }
// //   #endif
// //   return nc;
// // }
static inline coord normal_contact (coord ns, coord nf, double angle)
{
  coord nc; 
  #if dimension == 2
  if (- ns.x*nf.y + ns.y*nf.x > 0) { // 2D cross product // difference between pointing towards or away
    
    //  standard matrix *-1
    nc.x = - ns.x*cos(angle) + ns.y*sin(angle);
    nc.y = - ns.x*sin(angle) - ns.y*cos(angle);
  }
  else { // inverted  this is for axi-part
    nc.x = - ns.x*cos(angle) - ns.y*sin(angle);
    nc.y = + ns.x*sin(angle) - ns.y*cos(angle);
  }
 
  #elif dimension == 3
  if(1==0){
          coord rot_axis;
          // normal vector of ns and nf calculated by cross product.
          rot_axis.x = -ns.y*nf.z + ns.z*nf.y;
          rot_axis.y = -ns.z*nf.x + ns.x*nf.z;
          rot_axis.z = -ns.x*nf.y + ns.y*nf.x;
          normalize(&rot_axis);
          // rotate ns by angle around rot_axis
          nc.x = - ns.x*(sq(rot_axis.x)          *(1.-cos(angle)) +             cos(angle)) 
                - ns.y*(rot_axis.x*rot_axis.y*(1.-cos(angle)) - rot_axis.z*sin(angle)) 
                - ns.z*(rot_axis.x*rot_axis.z*(1.-cos(angle)) + rot_axis.y*sin(angle));
              
          nc.y = - ns.x*(rot_axis.y*rot_axis.x*(1.-cos(angle)) + rot_axis.z*sin(angle)) 
                - ns.y*(sq(rot_axis.y)       *(1.-cos(angle)) +            cos(angle)) 
                - ns.z*(rot_axis.y*rot_axis.z*(1.-cos(angle)) -rot_axis.x*sin(angle));
              
          nc.z = - ns.x*(rot_axis.z*rot_axis.x*(1.-cos(angle)) - rot_axis.y*sin(angle)) 
                - ns.y*(rot_axis.z*rot_axis.y*(1.-cos(angle)) + rot_axis.x*sin(angle)) 
                - ns.z*(sq(rot_axis.z)       *(1.-cos(angle)) +            cos(angle));  
  }else{
          coord rot_axis;
          coord mid_vector;
          double val1=0.0;;
          foreach_dimension(){
            val1=val1+nf.x*ns.x;
          }
          //mid_vecor nf - (nf*ns)*ns
          foreach_dimension(){
              mid_vector.x = nf.x - val1*ns.x;
          }
          normalize(&mid_vector);
          // normal vector of ns and nf calculated by cross product.
          // rot_axis.x = -ns.y*nf.z + ns.z*nf.y;
          // rot_axis.y = -ns.z*nf.x + ns.x*nf.z;
          // rot_axis.z = -ns.x*nf.y + ns.y*nf.x;
          rot_axis.x = -ns.y*mid_vector.z + ns.z*mid_vector.y;
          rot_axis.y = -ns.z*mid_vector.x + ns.x*mid_vector.z;
          rot_axis.z = -ns.x*mid_vector.y + ns.y*mid_vector.x;
          normalize(&rot_axis);
          // rotate ns by angle around rot_axis
          nc.x = - ns.x*(sq(rot_axis.x)          *(1.-cos(angle)) +             cos(angle)) 
                - ns.y*(rot_axis.x*rot_axis.y*(1.-cos(angle)) - rot_axis.z*sin(angle)) 
                - ns.z*(rot_axis.x*rot_axis.z*(1.-cos(angle)) + rot_axis.y*sin(angle));
              
          nc.y = - ns.x*(rot_axis.y*rot_axis.x*(1.-cos(angle)) + rot_axis.z*sin(angle)) 
                - ns.y*(sq(rot_axis.y)       *(1.-cos(angle)) +            cos(angle)) 
                - ns.z*(rot_axis.y*rot_axis.z*(1.-cos(angle)) -rot_axis.x*sin(angle));
              
          nc.z = - ns.x*(rot_axis.z*rot_axis.x*(1.-cos(angle)) - rot_axis.y*sin(angle)) 
                - ns.y*(rot_axis.z*rot_axis.y*(1.-cos(angle)) + rot_axis.x*sin(angle)) 
                - ns.z*(sq(rot_axis.z)       *(1.-cos(angle)) +            cos(angle));  
  }
  #endif
  return nc;
}

/**
This function is an adaptation of the
[reconstruction()](/src/fractions.h#reconstruction) function which
takes into account the contact angle on embedded boundaries. */

void reconstruction_contact (scalar f, vector n, scalar alpha, scalar intersect_true)
{

  /**
  We first reconstruct the (n, alpha) fields everywhere, using the
  standard function. */
  
scalar f_middle[];
 foreach(){
    f_middle[] = f[];
 }
 foreach(){
    if(cs[]<=0.0){
       f[] = f_height[];
    }
 }
 //reconstruction n,alpha from f, but correct set the ghost volume fraction, 
 //but we could not direct give ghost volume fraction to f[], sin it will change
 //volume of mixed cell.
 //f_height in reconstruction4 is only a temp scalar.. 
 //this newly build n and alpha, could not use to calculate fraction cs[]>0 for f.  

 if(1==0){
  reconstruction4 (f, n, alpha); //height function->normal interface_normal --- mycs (point, c)
 }else{
  reconstruction (f, n, alpha); //MYC approximation 
 }

 foreach(){
   f[] = f_middle[];
 }

  /**
  In cells which contain an embedded boundary and an interface, we
  modify the reconstruction to take the contact angle into account. */
  
  foreach(){
   // printf("ddddddd\n");
    if (cs[] < 1. && cs[] > 0. && f[] < 1.  && f[] > 0.) {
        // coord nsf = (coord){solid_n.x[],solid_n.y[]};
        // double alphasf = solid_neg_alpha[];
        // coord nlg = (coord){n.x[],n.y[]};
        // double alphalg=alpha[];
        // double a1,b1,c1,a2,b2,c2;
        // a1 = nsf.x,b1=nsf.y,c1=-0.5*(nsf.x+nsf.y)-alphasf;
        // a2 = nlg.x,b2=nlg.y,c2=-0.5*(nlg.x+nlg.y)-alphalg; 
        // bool intersect_flag= is_intersection_within_unit_square(a1,b1,c1,a2,b2,c2);
        // if(intersect_flag){
          if(intersect_true[]==1){
                coord ns = facet_normal (point, cs, fs);
                coord nf;
                foreach_dimension()
                    nf.x = n.x[];
              //  printf("ns: %g %g %g\n",ns.x,ns.y,ns.z);
                normalize (&ns);
                coord nc = normal_contact (ns, nf, contact_angle[]);
                foreach_dimension()
                  n.x[] = nc.x;

              //  printf("n: %g %g %g\n",n.x[],n.y[],n.z[]);
                alpha[] = line_alpha (f[], nc);
        }
    }
    boundary ({n, alpha});  
  }
}

typedef struct {
    int id;
    double x, y, alpha, n_x, n_y, c1,c2,sn_x,sn_y;
    double triplex,tripley;
} CellInfo;

/**
At every timestep, we modify the volume fraction values in the
embedded solid to take the contact angle into account. */

bool is_intersection_within_unit_square(double a1, double b1, double c1, double a2, double b2, double c2, double* xx, double* yy) {
    double det = a1 * b2 - a2 * b1;
    // Point intersection;

    // Check if lines are not parallel
    if (fabs(det) > 1e-9) {
        double intersect_x,intersect_y;
        intersect_x = (b1 * c2 - b2 * c1) / det;
        intersect_y = (a2 * c1 - a1 * c2) / det;

        // Check if intersection is within the unit square
        if (0 <= intersect_x && intersect_x <= 1 && 0 <= intersect_y && intersect_y <= 1) {
            *xx = intersect_x;// - 0.5;
            *yy = intersect_y;// - 0.5;
            return true;
        }
    }
    return false;
}




event contact (i++)
{
  vector n[];
  scalar alpha[];

  /**
  We first reconstruct (n,alpha) everywhere. Note that this is
  necessary since we will use neighborhood stencils whose values may
  be reconstructed by adaptive and/or parallel boundary conditions. */
  reconstruction4 (ff, n, alpha);
  scalar cs_neg[];
  foreach(){
    cs_neg[] = 1.0-cs[];
  }
  scalar solid_neg_alpha[];
  vector solid_n[];
  reconstruction (cs_neg, solid_n, solid_neg_alpha);
  // scalar intersect_true[];
  vector triple_coord[];
  int total_number_triple=0;

  scalar intersect_true_temp[];

  foreach(){
    intersect_true_temp[] = intersect_true[];
  }

  foreach(reduction(+:total_number_triple)){
        intersect_true[] = 0;
        f_height[] = ff[];
        foreach_dimension(){
          triple_coord.x[] = 0;

        }
      // printf("ddddddd\n");
        if (cs[] < 1. && cs[] > 0. && ff[] < 1.  && ff[] > 0.) {
            coord nsf = (coord){solid_n.x[],solid_n.y[]};
            double alphasf = solid_neg_alpha[];
            coord nlg = (coord){n.x[],n.y[]};
            double alphalg=alpha[];
            double a1,b1,c1,a2,b2,c2;
            a1 = nsf.x,b1=nsf.y,c1=-0.5*(nsf.x+nsf.y)-alphasf;
            a2 = nlg.x,b2=nlg.y,c2=-0.5*(nlg.x+nlg.y)-alphalg; 
            double xx=0,yy=0;
            bool intersect_flag= is_intersection_within_unit_square(a1,b1,c1,a2,b2,c2,&xx,&yy);
            if(intersect_flag){
                intersect_true[] = 1;
                total_number_triple+=1;
                triple_coord.x[] = (xx-0.5)*Delta+x;
                triple_coord.y[] = (yy-0.5)*Delta+y;
            }
        }
  }
  
  bool flag33=false;
  int method=0;
  if(total_number_triple==0 && 1==0){
      flag33=true;
      foreach(reduction(+:total_number_triple)){
         if(intersect_true_temp[]==1 && level==level_interface){
              if(ff[]>0.0 && ff[]<1.0 && cs[]<1.0 && cs[]>0.0 && intersect_true[]==0){
                  intersect_true[]=1;
                  total_number_triple+=1;
              }
         }
      }
      if(total_number_triple>=1){
            method=1;
      }
      if(total_number_triple==0){
            foreach(reduction(+:total_number_triple)){
                if(cs[]>0.0 && cs[]<1.0 && ff[]>0.0 && ff[]<1.0 && level==level_interface){
                  bool flag1 = false;
                  foreach_neighbor(1){
                      if(intersect_true_temp[]==1){
                          flag1 = true;
                      }
                  }
                  if(flag1 && intersect_true[]==0){
                     intersect_true[]=1;
                      total_number_triple+=1;
                  }
                }
            }
            if(total_number_triple>=1){
                  method=2;
            }
      }

      // if(total_number_triple==0){
      //       foreach(reduction(+:total_number_triple)){
      //           if(cs[]>0.0 && cs[]<1.0 && ff[]>0.0 && ff[]<1.0 && level==level_interface){
      //             bool flag1 = false;
      //             foreach_neighbor(1){
      //                 if(topo_mask_s[]==1 && ff[]>0.0 && ff[]<1.0){
      //                     flag1 = true;
      //                 }
      //             }
      //             if(flag1 && intersect_true[]==0){
      //                intersect_true[]=1;
      //                 total_number_triple+=1;
      //             }
      //           }
      //       }
      //       if(total_number_triple>=1){
      //               method=3;
      //         }
      // }
  }
 if(total_number_triple<=0){
                    method=-1;
  }

bool flag34=false;
if(!flag33 && 1==0){
  flag34=true;
  scalar intersect_true_temp2[];
  foreach(reduction(+:total_number_triple)){
    intersect_true_temp2[]=0;
     if(level==level_interface && intersect_true_temp[]==1 && intersect_true[]==0){
        bool flag=false;
        foreach_neighbor(1){
           if(intersect_true[]==1){
              flag=true;
           }
        }
        if((!flag) && intersect_true[]==0){
            intersect_true_temp2[] =2;
            // total_number_triple+=1;
        }
     }
  }
  foreach(reduction(+:total_number_triple)){
      if(level==level_interface && cs[]>0.0 && cs[]<1.0 && ff[]>0.0 && ff[]<1.0 && intersect_true[]==0){
         bool flag=false;
         Point me;
         foreach_neighbor(1){
            if(!(me.i==point.i && me.j==point.j)){
              if(intersect_true_temp2[]==2){
                  flag=true;
              }
            }
         }
         if(flag && intersect_true[]==0){
            intersect_true[]=2;
            total_number_triple+=1;
         }
      }
  }
}

if(1==1){
    foreach(reduction(+:total_number_triple)){
      if(intersect_true[]==1 || intersect_true[]==2){
          bool flag=false;
          double value = css_test[];
          foreach_neighbor(2){
            if(intersect_true[]==1 || intersect_true[]==2){
              if(value<css_test[]){
                flag=true;
              }
            }
          }
           if(flag){
                intersect_true[] = 0;
                total_number_triple -=1;
           }
      }
    }
}else{
    scalar mindist[];
    foreach(){
      if(intersect_true[]==1){
          bool flag=false;
          double mindist_local=HUGE;
          coord current;
          current.x = x, current.y=y;
          foreach_neighbor(2){
            coord here;
            if(topo_mask_s[]==0 && ff[]<=0){
              here.x = x, here.y = y;
              double distance = sqrt(sq(current.x-here.x) + sq(current.y-here.y));
              if(distance < mindist_local){
                    mindist_local = distance;
              }
            }
          }
          mindist[] = mindist_local;
      }
    }

   foreach(reduction(+:total_number_triple)) {
        if(intersect_true[]==1){
            bool flag=false;
            double value = mindist[];
            foreach_neighbor(2){
              if(intersect_true[]==1){
                if(value<mindist[]){
                  flag=true;
                }
              }
            }
            if(flag){
                intersect_true[] = 0;
                total_number_triple -=1;
            }
        }
    }
}


  
  
 


  if(flag33 && total_number_triple>=1){ //calculate triple using (-nsy,nsx)
      foreach(){
      // printf("ddddddd\n");
        if (intersect_true[]==1) {
            coord nsf = (coord){solid_n.x[],solid_n.y[]};
            double alphasf = solid_neg_alpha[];
            // coord nlg = (coord){n.x[],n.y[]};
            coord nlg = (coord){-solid_n.y[],solid_n.x[]}; //using (-ny,nx) ->must have intersection
            bool same_direction= (n.x[]*(-solid_n.y[])+n.y[]*solid_n.x[])>0;
            if(!same_direction){
                foreach_dimension(){
                   nlg.x = - nlg.x;
                }
            }
            // double alphalg=alpha[];
            double alphalg =plane_alpha (ff[], nlg);//recalculate alphalg            
            double a1,b1,c1,a2,b2,c2;
            a1 = nsf.x,b1=nsf.y,c1=-0.5*(nsf.x+nsf.y)-alphasf;
            a2 = nlg.x,b2=nlg.y,c2=-0.5*(nlg.x+nlg.y)-alphalg; 
            double xx=0,yy=0;
            bool intersect_flag= is_intersection_within_unit_square(a1,b1,c1,a2,b2,c2,&xx,&yy);
            if(intersect_flag){
                intersect_true[] = 1;
                triple_coord.x[] = (xx-0.5)*Delta+x;
                triple_coord.y[] = (yy-0.5)*Delta+y;
            }else{
                intersect_true[] = 0;
            }
        }
      }
      
  }else if(flag34 && total_number_triple>=1){
      foreach(){
      // printf("ddddddd\n");
        if (intersect_true[]==2) {
            coord nsf = (coord){solid_n.x[],solid_n.y[]};
            double alphasf = solid_neg_alpha[];
            // coord nlg = (coord){n.x[],n.y[]};
            coord nlg = (coord){-solid_n.y[],solid_n.x[]}; //using (-ny,nx) ->must have intersection
            bool same_direction= (n.x[]*(-solid_n.y[])+n.y[]*solid_n.x[])>0;
            if(!same_direction){
                foreach_dimension(){
                   nlg.x = - nlg.x;
                }
            }
            // double alphalg=alpha[];
            double alphalg =plane_alpha (ff[], nlg);//recalculate alphalg            
            double a1,b1,c1,a2,b2,c2;
            a1 = nsf.x,b1=nsf.y,c1=-0.5*(nsf.x+nsf.y)-alphasf;
            a2 = nlg.x,b2=nlg.y,c2=-0.5*(nlg.x+nlg.y)-alphalg; 
            double xx=0,yy=0;
            bool intersect_flag= is_intersection_within_unit_square(a1,b1,c1,a2,b2,c2,&xx,&yy);
            if(intersect_flag){
                intersect_true[] = 1;
                triple_coord.x[] = (xx-0.5)*Delta+x;
                triple_coord.y[] = (yy-0.5)*Delta+y;
            }else{
                intersect_true[] = 0;
            }
        }
      }
  }
  

reconstruction_contact (ff, n, alpha, intersect_true);
///////////////for larger region
int number = 0;
foreach (reduction(+:number)) {
    if (intersect_true[] == 1) {
        number++;
        // printf("number %d :x=%g, y=%g process=%d\n",number, x,y,pid());
    }
}

if(number==0 && pid()==0){
    printf("number of triple point ==0\n");
}

int selectedCount = 0;
// scalar tag_cell[];
foreach (serial) {
    // tag_cell[] = 0;
    if (intersect_true[] == 1) {
        if (selectedCount < number) {  // Ensure we don't exceed the array size
            selectedCount++;
        }
    }
}
int number_local = selectedCount;
// CellInfo selectedCells[number_local];
CellInfo* selectedCells = (CellInfo*)malloc(number_local * sizeof(CellInfo));
int globalCount;
MPI_Allreduce(&selectedCount, &globalCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
int offset;
MPI_Scan(&selectedCount, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
offset -= selectedCount;  // Subtract local count to get the offset
int local_number = 0;
foreach (serial) {
    if (intersect_true[] == 1) {
        if (local_number < selectedCount) {
            int globalID = offset + local_number; // Corrected: Use local_number here
            selectedCells[local_number].id = globalID;
            // tag_cell[] = globalID;
            selectedCells[local_number].x = x;//globalID * 0.1;
            selectedCells[local_number].y = y;//globalID * 0.2;
            selectedCells[local_number].alpha = alpha[];//globalID * 0.3;
            selectedCells[local_number].n_x = n.x[];
            selectedCells[local_number].n_y = n.y[];//globalID * 0.5;
            selectedCells[local_number].c1 = cs[];
            selectedCells[local_number].c2 = ff[];
            selectedCells[local_number].sn_x = solid_n.x[];
            selectedCells[local_number].sn_y = solid_n.y[];//globalID * 0.5;
            selectedCells[local_number].triplex = triple_coord.x[];
            selectedCells[local_number].tripley = triple_coord.y[];//globalID * 0.5;
            local_number++;
            // printf("x=%g,y=%g\n",x,y);
        }
    }
}
// Step 1: Gather counts to root process
        int *allCounts = NULL;
        if (pid() == 0) {
            allCounts = (int*)malloc(npe() * sizeof(int));
        }
        MPI_Gather(&selectedCount, 1, MPI_INT, allCounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Step 2: Root process allocates space for all selected cells
        CellInfo *allSelectedCells = NULL;
        int totalSelectedCount = 0;
        int *displs = NULL;
        if (pid() == 0) {
            for (int i = 0; i < npe(); i++) {
                totalSelectedCount += allCounts[i];
            }
            allSelectedCells = (CellInfo*)malloc(totalSelectedCount * sizeof(CellInfo));
            displs = (int*)malloc(npe() * sizeof(int));
            displs[0] = 0;
            for (int i = 1; i < npe(); i++) {
                displs[i] = displs[i - 1] + allCounts[i - 1];
            }
        }

        // Step 3: Create a custom MPI Datatype for the structure
        MPI_Datatype cellInfoType;
        // MPI_Type_contiguous(6, MPI_DOUBLE, &cellInfoType);
        // MPI_Type_contiguous(8, MPI_DOUBLE, &cellInfoType);
        MPI_Type_contiguous(12, MPI_DOUBLE, &cellInfoType);
        MPI_Type_commit(&cellInfoType);

        // Step 4: Gather selected cells at root process
        MPI_Gatherv(selectedCells, selectedCount, cellInfoType,
                    allSelectedCells, allCounts, displs, cellInfoType,
                    0, MPI_COMM_WORLD);

        // Step 5: Broadcast the complete array to all processes
        MPI_Bcast(&totalSelectedCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (pid() != 0) {
            allSelectedCells = (CellInfo*)malloc(totalSelectedCount * sizeof(CellInfo));
        }
        MPI_Bcast(allSelectedCells, totalSelectedCount, cellInfoType, 0, MPI_COMM_WORLD);

        // Example: Each process prints its received data (for demonstration)
// printf("second-22222222222222\n");
        // printf("number=%d\n",number);
        // printf("Process %d received data:\n", pid());
        // for (int i = 0; i < totalSelectedCount; ++i) {
        //     printf("Cell %d: x = %f, y = %f\n", allSelectedCells[i].id, allSelectedCells[i].x, allSelectedCells[i].y);
        // }



  /**
  We then look for "contact" cells in the neighborhood of each cell
  entirely contained in the embedded solid. */

scalar corner_ff_temp[];
if(1==0){
    foreach(){
        corner_ff_temp[]=0;
        corner_ff[]=0;
        if(cs[]<1.0 && cs[]>0.0 && level==level_interface){
            int number=0;
            Point me = point;
            foreach_neighbor(1){
                if(!(me.i==point.i && me.j==point.j)){
                    if(cs[]>=1){
                        number+=1;
                    }
                }
            }
            if(number>=5){
                corner_ff_temp[] =1;
            }
        }
            
    }
    foreach(){
      corner_ff[]=0;
      if(cs[]<1.0 && cs[]>0.0){
          bool flag=false;
          foreach_neighbor(2){
                if(corner_ff_temp[]==1){
                    flag=true;
                }
          }
          if(flag){
              corner_ff[] = 1;
          }
      }
    }
    foreach(){
      corner_ff_temp[]=corner_ff[];
    }
    foreach(){
      if(cs[]<1.0 && cs[]>0.0 && corner_ff[]==0){
          bool flag=false;
          foreach_neighbor(2){
                if(corner_ff_temp[]==1){
                    flag=true;
                }
          }
          if(flag){
              corner_ff[] = 1;
          }
      }
    }
}else{
    foreach(){
      if(level==level_interface && topo_mask_s[]==0){
          corner_ff[]=1;
      }
    }
}


foreach(){
  //  region[]=0;
  f_height[] = ff[];
  bool flag=false;
  flag = (cs[]>0.0 && cs[]<1e-6) || ((cs[]<1.0) && (cs[]>1.0-1e-6));
  //  if (cs[] == 0. && (level==level_interface)) { //level>=level_interface
   if (cs[]<1.0 && (!flag) && (level==level_interface) && corner_ff[]==0) { 
      double fc = 0., sfc = 0.;
      coord o = {x, y, z};
      for(int ii=0;ii<totalSelectedCount;ii++){
        double fc_temp=0.0;//sfc_temp=0.0;
        bool insolid_flag=false;
        coord solidn;
          solidn.x = allSelectedCells[ii].sn_x;
          solidn.y = allSelectedCells[ii].sn_y;
          coord a;
          a.x = allSelectedCells[ii].x;
          a.y = allSelectedCells[ii].y;
          // printf("inside foreach arrayx=%g,arrayy=%g\n",a.x,a.y);
          double distance1 = fabs(a.x-x)/Delta;
          double distance2 = fabs(a.y-y)/Delta;
          if(distance1<0.2 && distance2<0.2){
                   // if(1==1 && intersect_true[] == 1){
                    if(1==0 && intersect_true[] == 1){
                      coord nf;
                      nf.x = allSelectedCells[ii].n_x;
                      nf.y = allSelectedCells[ii].n_y;
                      double alpha1 = allSelectedCells[ii].alpha;
                      coord triple_location;
                      triple_location.x =  allSelectedCells[ii].triplex;
                      triple_location.y =  allSelectedCells[ii].tripley;
                                  //modify alpha1 got from n and triple point


                                  //interface at triple point using nf_scale and alpha_scale
                      double nf_abs=0.0;// = fabs(nf.x) + fabs(nf.y) + fabs(nf.z);
                      foreach_dimension(){
                        nf_abs+=fabs(nf.x);
                      }
                      coord nf_scale;
                      foreach_dimension(){
                          nf_scale.x = nf.x/nf_abs;
                      }
                      double alpha_scale=0.0;
                      foreach_dimension(){
                          alpha_scale += (nf_scale.x*(triple_location.x-a.x)/Delta);
                      }
                      if(1==1){
                        foreach_dimension(){
                          nf.x = nf_scale.x;
                        }
                          alpha1 = alpha_scale;
                      }
                      f_height[] = line_area(nf.x, nf.y, alpha1);
                  }
           }else if(distance1<=4 && distance2<=4){
                  // region[] = 1;
                  coord direction;
                  coord triple_location;
                  triple_location.x =  allSelectedCells[ii].triplex;
                  triple_location.y =  allSelectedCells[ii].tripley;
                  double c1 = allSelectedCells[ii].c1; //cs
                  double c2 = allSelectedCells[ii].c2; //f
                  double coef = c1*(1. - c1)*c2*(1. - c2);
                  // sfc += coef;

                  coord nf;
                    nf.x = allSelectedCells[ii].n_x;
                    nf.y = allSelectedCells[ii].n_y;
                  double alpha1 = allSelectedCells[ii].alpha;
                  double nf_abs=0.0;// = fabs(nf.x) + fabs(nf.y) + fabs(nf.z);
                    foreach_dimension(){
                      nf_abs+=fabs(nf.x);
                    }
                    coord nf_scale;
                    foreach_dimension(){
                      nf_scale.x = nf.x/nf_abs;
                    }
                    double alpha_scale=0.0;
                    foreach_dimension(){
                      alpha_scale += (nf_scale.x*(triple_location.x-a.x)/Delta);
                    }
                  if(1==1){
                      foreach_dimension(){
                        nf.x = nf_scale.x;
                      }
                      alpha1 = alpha_scale;
                  }
                  coord b;
                  foreach_dimension()
                    a.x = (o.x - a.x)/Delta - 0.5, b.x = a.x + 1.;
                  // fc += coef*rectangle_fraction (nf, alpha1, a, b);
                  fc_temp = rectangle_fraction (nf, alpha1, a, b);  // extrapolated from center of the cell
                  double alpha_temp;
                  alpha_temp = line_alpha (fc_temp, nf);
                  coord p_temp;
                  double area1 = line_length_center (nf, alpha_temp, &p_temp);
                  p_temp.x = x + Delta*p_temp.x;
                  p_temp.y = y + Delta*p_temp.y;
                  direction.x = p_temp.x - triple_location.x;
                  direction.y = p_temp.y - triple_location.y;
                  insolid_flag = ((direction.x*solidn.x+direction.y*solidn.y)<0.0 && !(fabs(area1)<1e-8)); //1e-6 ;
                  if(cs[]==0){
                        insolid_flag = true;
                    }
                    if(insolid_flag){  
                          fc += coef*fc_temp;
                          sfc += coef;
                    }
                  
          }
      }
      if (sfc > 0.){
	      f_height[] = fc/sfc;
        if(cs[]==0){
            ff[] = fc/sfc;
        }
      }
   }
}
 boundary ({ff});

      free(selectedCells);
        // Finalize: Free allocated memory and finalize MPI
        if (pid() == 0) {
            free(allCounts);
            free(displs);
        }
        free(allSelectedCells);
        MPI_Type_free(&cellInfoType);
}

// fixme: adaptive mesh has to be done
