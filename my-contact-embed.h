/**
# Contact angles on an embedded boundary

This header file implemented contact angles for [VOF
interfaces](/src/vof.h) on [embedded boundaries](/src/embed.h).

The contact angle is defined by this field and can be constant or
variable. */

#include "./reconstruction4.h" //reconstruct using height_normal

extern scalar intersect_true;

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
  
 // reconstruction (f, n, alpha); // interface_normal --- mycs (point, c)
  reconstruction4 (f, n, alpha); // interface_normal --- mycs (point, c)

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

/**
At every timestep, we modify the volume fraction values in the
embedded solid to take the contact angle into account. */

bool is_intersection_within_unit_square(double a1, double b1, double c1, double a2, double b2, double c2) {
    double det = a1 * b2 - a2 * b1;
    Point intersection;

    // Check if lines are not parallel
    if (fabs(det) > 1e-9) {
        double intersect_x,intersect_y;
        intersect_x = (b1 * c2 - b2 * c1) / det;
        intersect_y = (a2 * c1 - a1 * c2) / det;

        // Check if intersection is within the unit square
        if (0 <= intersect_x && intersect_x <= 1 && 0 <= intersect_y && intersect_y <= 1) {
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
  foreach(){
        intersect_true[] = 0;
      // printf("ddddddd\n");
        if (cs[] < 1. && cs[] > 0. && ff[] < 1.  && ff[] > 0.) {
            coord nsf = (coord){solid_n.x[],solid_n.y[]};
            double alphasf = solid_neg_alpha[];
            coord nlg = (coord){n.x[],n.y[]};
            double alphalg=alpha[];
            double a1,b1,c1,a2,b2,c2;
            a1 = nsf.x,b1=nsf.y,c1=-0.5*(nsf.x+nsf.y)-alphasf;
            a2 = nlg.x,b2=nlg.y,c2=-0.5*(nlg.x+nlg.y)-alphalg; 
            bool intersect_flag= is_intersection_within_unit_square(a1,b1,c1,a2,b2,c2);
            if(intersect_flag){
                intersect_true[] = 1;
            }
        }
  }
  
  reconstruction_contact (ff, n, alpha, intersect_true);

  /**
  We then look for "contact" cells in the neighborhood of each cell
  entirely contained in the embedded solid. */
  
  foreach() {
    if (cs[] == 0.) {
      double fc = 0., sfc = 0.;
      coord o = {x, y, z};
      foreach_neighbor()
	if (cs[] < 1. && cs[] > 0. && ff[] < 1.  && ff[] > 0.) {	  
      if(intersect_true[]==1){
            /**
            This is a contact cell. We compute a coefficient meant to
            weight the estimated volume fraction according to how well
            the contact point is defined. We assume here that a contact
            point is better reconstructed if the values of `cs` and `f`
            are both close to 0.5. */

            double coef = cs[]*(1. - cs[])*ff[]*(1. - ff[]);
            sfc += coef;

            /**
            We then compute the volume fraction of the solid cell
            (centered on `o`), using the extrapolation of the interface
            reconstructed in the contact cell. */

            coord nf;
            foreach_dimension()
              nf.x = n.x[];
            coord a = {x, y, z}, b;
            foreach_dimension()
              a.x = (o.x - a.x)/Delta - 0.5, b.x = a.x + 1.;
            fc += coef*rectangle_fraction (nf, alpha[], a, b);
      }
	}

      /**
      The new volume fraction value of the solid cell is the weighted
      average of the volume fractions reconstructed from all
      neighboring contact cells. */
      
      if (sfc > 0.)
	ff[] = fc/sfc;
    }
  }
  boundary ({ff});
}

// fixme: adaptive mesh has to be done
