extern scalar css_test3,css_test3_n;
extern face vector fss_test3,fss_test3_n;
extern scalar ff;
// extern double Tks;
extern scalar css_test,css_test2;
extern face vector fss_test,fss_test2;
extern double Tsat00,Tsat0;

#include "fractions.h"
#include "curvature.h"

//extern face vector D;
// // //extern face vector cond_k;

event init(t=0){
    
}

// // // //#define is_solid(vv) vv>0.0 
// // // //if not specified alpha in multi-grid is  restriction_average
// // // //version1;
// // // void get_D(face vector D){
// // //     //  foreach_face(){
// // //     //     //fixed me: fss_test3 may be better, same situation in two-phase.h
// // //     //     if((!(is_solid(css_test3[]))) && (!(is_solid(css_test3[1]))) ){ //gas and fliud field
// // //     //         D.x[] = cond_k.x[];
// // //     //     }else if((is_solid(css_test3[])) && (is_solid(css_test3[1]))){ // face between solid
// // //     //         D.x[] = Tks;
// // //     //     }else{ //face between solid and fluid; cond_k = k_mixed;
// // //     //         D.x[] = 2.0*cond_k.x[]*Tks/(cond_k.x[]+Tks); // only_for_solid_extrictly_lie_on_the_face_of_cell
// // //     //     }
// // //     //  }
// // //     foreach_face(face vector D){
// // //          D.x[] = cond_k.x[];
// // //     }
// // // }

// // // //version2, solid and fluid could be in a cell.
// // // void get_D2(face vector D){
// // //      foreach_face(){
// // //         //fixed me: fss_test3 may be better, same situation in two-phase.h
// // //         if(fss_test3_n.x[]>=1.0){ //gas and fliud field
// // //             D.x[] = cond_k.x[];
// // //         }else if(fss_test3_n.x[]<=0.0){ // face between solid
// // //             D.x[] = Tks;
// // //         }else{ //fixed be: separately caculate fluid and solid not developed.
// // //             //printf("error in computing D2, in present version solid occupy the whole cell, so fss_test3_n=0or1\n");
// // //             printf("need to developed,fluid and solid calculate sepatately");
// // //             exit(1);
// // //             //D.x[] = 2.0*cond_k.x[]*Tks/(cond_k.x[]+Tks)*1.0; // only_for_solid_extrictly_lie_on_the_face_of_cell
// // //         }
// // //      }
// // // }

void embed_fraction_refine_s (Point point, scalar css_test3)
{
  double cc = css_test3[];

  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (cc <= 0. || cc >= 1.) {
    foreach_child()
      css_test3[] = cc;
  }
  else {

    /**
    If the cell contains the embedded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, css_test3, fss_test3);
    double alpha = plane_alpha (cc, n);
      
    foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      foreach_dimension()
	nc.x = child.x*n.x;
      css_test3[] = rectangle_fraction (nc, alpha, a, b);
    }
  }
}


void embed_fraction_refine_s_n (Point point, scalar css_test3_n)
{
  double cc = css_test3_n[];

  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (cc <= 0. || cc >= 1.) {
    foreach_child()
      css_test3_n[] = cc;
  }
  else {

    /**
    If the cell contains the embedded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, css_test3_n, fss_test3_n);
    double alpha = plane_alpha (cc, n);
      
    foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      foreach_dimension()
	nc.x = child.x*n.x;
      css_test3_n[] = rectangle_fraction (nc, alpha, a, b);
    }
  }
}


//refine for fss_test3

foreach_dimension()
void embed_face_fraction_refine_s_x (Point point, scalar s)
{
  vector fss_test3 = s.v;

  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (css_test3[] <= 0. || css_test3[] >= 1.) {

    /**
    We need to make sure that the fine cells face fractions match
    those of their neighbours. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	fine(fss_test3.x,1,j,k) = css_test3[];
    for (int i = 0; i <= 1; i++)
      if (!is_refined(neighbor(2*i-1)) && neighbor(2*i-1).neighbors &&
	  (is_local(cell) || is_local(neighbor(2*i-1))))
	for (int j = 0; j <= 1; j++)
	  for (int k = 0; k <= 1; k++)
	    fine(fss_test3.x,2*i,j,k) = fss_test3.x[i];
  }
  else {

    /**
    If the cell contains the embedded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, css_test3, fss_test3);
    double alpha = plane_alpha (css_test3[], n);
      
    /**
    We need to reconstruct the face fractions *fss_test3* for the fine cells.
    
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
      fine(fss_test3.x,1,1 - i) = n.y < 0. ? 1. - i : i;
      fine(fss_test3.x,1,i) = n.y < 0. ? i - 2.*yc : 1. - i + 2.*yc;
    }
    else
      fine(fss_test3.x,1,0) = fine(fss_test3.x,1,1) = alpha > 0.;

#else // dimension == 3

    /**
    in 3D, we use the 2D projection of the reconstruction. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	if (!fine(css_test3,0,j,k) || !fine(css_test3,1,j,k))
	  fine(fss_test3.x,1,j,k) = 0.;
	else {
	  static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
	  coord nc;
	  nc.x = 0., nc.y = (2.*j - 1.)*n.y, nc.z = (2.*k - 1.)*n.z;
	  fine(fss_test3.x,1,j,k) = rectangle_fraction (nc, alpha, a, b);
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
	  if (fss_test3.x[i] <= 0. || fss_test3.x[i] >= 1.)
	    for (int j = 0; j <= 1; j++)
	      for (int k = 0; k <= 1; k++)
		fine(fss_test3.x,2*i,j,k) = fss_test3.x[i];
	  else {
#if dimension == 2
	  
	    /**
	    In 2D the orientation is obtained by looking at the values
	    of face fractions in the transverse direction. */
	  
	    double a = fss_test3.y[0,1] <= 0. || fss_test3.y[2*i-1,1] <= 0. ||
	      fss_test3.y[] >= 1. || fss_test3.y[2*i-1] >= 1.;
	    if ((2.*a - 1)*(fss_test3.x[i] - 0.5) > 0.) {
	      fine(fss_test3.x,2*i,0) = a;
	      fine(fss_test3.x,2*i,1) = 2.*fss_test3.x[i] - a;
	    }
	    else {
	      fine(fss_test3.x,2*i,0) = 2.*fss_test3.x[i] + a - 1.;
	      fine(fss_test3.x,2*i,1) = 1. - a;
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
		fine(fss_test3.x,2*i,j,k) =
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
	    if (fine(fss_test3.x,2*i,j,k) && !fine(css_test3,i,j,k))
	      fine(fss_test3.x,2*i,j,k) = 0.;
      }
  }
}



//refine for fss_test3_n


// foreach_dimension()
// void embed_face_fraction_refine_s_n_x (Point point, scalar s)
foreach_dimension()
void embed_face_fraction_refine_s_n_x (Point point, scalar s)
{
  vector fss_test3_n = s.v;

  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (css_test3_n[] <= 0. || css_test3_n[] >= 1.) {

    /**
    We need to make sure that the fine cells face fractions match
    those of their neighbours. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	fine(fss_test3_n.x,1,j,k) = css_test3_n[];
    for (int i = 0; i <= 1; i++)
      if (!is_refined(neighbor(2*i-1)) && neighbor(2*i-1).neighbors &&
	  (is_local(cell) || is_local(neighbor(2*i-1))))
	for (int j = 0; j <= 1; j++)
	  for (int k = 0; k <= 1; k++)
	    fine(fss_test3_n.x,2*i,j,k) = fss_test3_n.x[i];
  }
  else {

    /**
    If the cell contains the embedded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, css_test3_n, fss_test3_n);
    double alpha = plane_alpha (css_test3_n[], n);
      
    /**
    We need to reconstruct the face fractions *fss_test3_n* for the fine cells.
    
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
      fine(fss_test3_n.x,1,1 - i) = n.y < 0. ? 1. - i : i;
      fine(fss_test3_n.x,1,i) = n.y < 0. ? i - 2.*yc : 1. - i + 2.*yc;
    }
    else
      fine(fss_test3_n.x,1,0) = fine(fss_test3_n.x,1,1) = alpha > 0.;

#else // dimension == 3

    /**
    in 3D, we use the 2D projection of the reconstruction. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	if (!fine(css_test3_n,0,j,k) || !fine(css_test3_n,1,j,k))
	  fine(fss_test3_n.x,1,j,k) = 0.;
	else {
	  static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
	  coord nc;
	  nc.x = 0., nc.y = (2.*j - 1.)*n.y, nc.z = (2.*k - 1.)*n.z;
	  fine(fss_test3_n.x,1,j,k) = rectangle_fraction (nc, alpha, a, b);
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
	  if (fss_test3_n.x[i] <= 0. || fss_test3_n.x[i] >= 1.)
	    for (int j = 0; j <= 1; j++)
	      for (int k = 0; k <= 1; k++)
		fine(fss_test3_n.x,2*i,j,k) = fss_test3_n.x[i];
	  else {
#if dimension == 2
	  
	    /**
	    In 2D the orientation is obtained by looking at the values
	    of face fractions in the transverse direction. */
	  
	    double a = fss_test3_n.y[0,1] <= 0. || fss_test3_n.y[2*i-1,1] <= 0. ||
	      fss_test3_n.y[] >= 1. || fss_test3_n.y[2*i-1] >= 1.;
	    if ((2.*a - 1)*(fss_test3_n.x[i] - 0.5) > 0.) {
	      fine(fss_test3_n.x,2*i,0) = a;
	      fine(fss_test3_n.x,2*i,1) = 2.*fss_test3_n.x[i] - a;
	    }
	    else {
	      fine(fss_test3_n.x,2*i,0) = 2.*fss_test3_n.x[i] + a - 1.;
	      fine(fss_test3_n.x,2*i,1) = 1. - a;
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
		fine(fss_test3_n.x,2*i,j,k) =
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
	    if (fine(fss_test3_n.x,2*i,j,k) && !fine(css_test3_n,i,j,k))
	      fine(fss_test3_n.x,2*i,j,k) = 0.;
      }
  }
}


//restriction_embed_linear;  
//restriction by T, not energy, parce que: gas temperature will not hugely affected by water energy


//without any special implimentation
void restriction_three_phase1 (Point point, scalar s)  //the same with restriction_volume_average,except cm
{
  double sum = 0.;
  foreach_child(){
    sum += s[];
    //sum += cm[]*s[];
  }
  s[] = sum/(1 << dimension);///(cm[] + 1e-30);

  //   double sum = 0.;
  // foreach_child()
  //   sum += cm[]*s[];
  // s[] = sum/(1 << dimension)/(cm[] + 1e-30);
}


double bilinear3 (Point point, scalar s)
//double bilinear3 (Point point, scalar s)
{
  #if dimension == 1
    return (3.*coarse(s) + coarse(s,child.x))/4.;
  #elif dimension == 2
    return (9.*coarse(s) + 
	    3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
	    coarse(s,child.x,child.y))/16.;
  #else // dimension == 3
    return (27.*coarse(s) + 
	    9.*(coarse(s,child.x) + coarse(s,0,child.y) +
		coarse(s,0,0,child.z)) + 
	    3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
		coarse(s,0,child.y,child.z)) + 
	    coarse(s,child.x,child.y,child.z))/64.;
  #endif
}

void refine_three_phase11 (Point point, scalar s){ //same with bilinear
   foreach_child(){
        s[] = bilinear3 (point, s);
   }
}

void refine_three_phase1 (Point point, scalar s){ //same with bilinear
   foreach_child(){
        #if dimension == 1
            s[] = (3.*coarse(s) + coarse(s,child.x))/4.;
        #elif dimension == 2
            s[] =  (9.*coarse(s) + 
                3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
                coarse(s,child.x,child.y))/16.;
        #else // dimension == 3
            s[] = (27.*coarse(s) + 
                9.*(coarse(s,child.x) + coarse(s,0,child.y) +
                coarse(s,0,0,child.z)) + 
                3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
                coarse(s,0,child.y,child.z)) + 
                coarse(s,child.x,child.y,child.z))/64.;
        #endif
   }
}


// attribute {
//   bool restriction_bool;
// }

// void restriction_three_phase2 (Point point, scalar s)  //the same with restriction_volume_average,except cm
// {
//   bool flag = s.restriction_bool;
//   double sum = 0.;
//   foreach_child(){
//     sum += s[];
//     //sum += cm[]*s[];
//   }
//   s[] = sum/(1 << dimension);///(cm[] + 1e-30);
// }

// void refine_three_phase2 (Point point, scalar s){ //same with bilinear
//    foreach_child(){
//         #if dimension == 1
//             s[] = (3.*coarse(s) + coarse(s,child.x))/4.;
//         #elif dimension == 2
//             s[] =  (9.*coarse(s) + 
//                 3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
//                 coarse(s,child.x,child.y))/16.;
//         #else // dimension == 3
//             s[] = (27.*coarse(s) + 
//                 9.*(coarse(s,child.x) + coarse(s,0,child.y) +
//                 coarse(s,0,0,child.z)) + 
//                 3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
//                 coarse(s,0,child.y,child.z)) + 
//                 coarse(s,child.x,child.y,child.z))/64.;
//         #endif
//    }
// }

// static inline void refine_bilinear (Point point, scalar s)
// {
//   foreach_child()
//     s[] = bilinear (point, s);
// }

// void restriction_three_phase1 (Point point, scalar s)
// {

//   // by now, I am not sure how to deal with one cell contain gas fliud and solid. 
//   if(!css_test3[] && css_test3_n[]){ //fluid and gas side
//             double val = 0., nv = 0.;
//             for (int i = 0; i <= 1; i++)
//             #if dimension > 2
//                 for (int j = 0; j <= 1; j++)
//             #endif
//                 if (fine(css_test3_n,0,i,j) && fine(css_test3_n,1,!i,!j))
//                 val += (fine(s,0,i,j) + fine(s,1,!i,!j))/2., nv++;
//             if (nv > 0.) {
//                 s[] = val/nv;
//                 return;
//             }

//             /**
//              Otherwise, we use the average of the child cells which are defined
//             (there is at least one). */
            
//             coord p = {0.,0.,0.};
//             foreach_child()
//                 if (css_test3_n[])
//                 p.x += x, p.y += y, p.z += z, val += s[], nv++;
//             assert (nv > 0.);
//             s[] = val/nv;
//   }else if(css_test3[] && !css_test3_n[]){  // for solid

//             /**
//              We first try to interpolate "diagonally". If enough child cells are
//             defined (i.e. have non-zero embedded fractions), we return the
//             corresponding value. */

//             double val = 0., nv = 0.;
//             for (int i = 0; i <= 1; i++)
//             #if dimension > 2
//                 for (int j = 0; j <= 1; j++)
//             #endif
//                 if (fine(css_test3,0,i,j) && fine(css_test3,1,!i,!j))
//                 val += (fine(s,0,i,j) + fine(s,1,!i,!j))/2., nv++;
//             if (nv > 0.) {
//                 s[] = val/nv;
//                 return;
//             }

//             /**
//              Otherwise, we use the average of the child cells which are defined
//             (there is at least one). */
            
//             coord p = {0.,0.,0.};
//             foreach_child()
//                 if (css_test3[]){
//                   p.x += x, p.y += y, p.z += z, val += s[], nv++;
//                 }
//             assert (nv > 0.);
//             s[] = val/nv;

//             // coord p = {0.,0.,0.};
//             // double weight=0.0;
//             // foreach_child()
//             //     if (css_test3[]){
//             //       val += s[]*rhocp[]*css_test3[], nv++;
//             //       weight += rhocp[]*css_test3[];
//             //     }
//             // assert (nv > 0.);
//             // s[] = val/weight;

//             // if (s.embed_gradient && s.boundary[0] != s.boundary_homogeneous[0]) {
//             //     coord o = {x,y,z}, g;
//             //     s.embed_gradient (point, s, &g);
//             //     foreach_dimension()
//             //     s[] += (o.x - p.x/nv)*g.x;
//             // }
//   }else if(css_test3[] && css_test3_n[]){
//       if(css_test3[]_n >= 0.5){ // this is fluid cell
//              //

//       }else{ // this is solid cell 


//       }






//       printf("within one cell: there are mixed(gas liquid) and solid, how to deal with this situation???????");
//       printf("error css_test3 and css_test3_n not consistent: css_test3=%g,css_test3_n=%g\n",css_test3[],css_test3_n[]);
//       exit(1);
//   }
// }

attribute {
  scalar f1; //fss_test3.x
  scalar f2; //fss_test.x
}
void D_restriction(Point point, scalar s){
    scalar f1 = s.f1;
    scalar f2 = s.f2;
}

// event init(t=0){
//       T.restriction = restriction_three_phase1;
//       T.refine = refine_three_phase1; // this should be consistent with function used in posson3.h; 
//       T.prolongation = refine_three_phase1;
//     //  T.dirty = true;
// }


//suppose we have already refine and restriction value
// foreach_dimension()
// void D_refine_x (Point point, scalar s){
//      vector v = s.v;
//      for (int j = 0; j <= 1; j++)
//       for (int k = 0; k <= 1; k++)
// 	fine(v.x,1,j,k) = 1;//css_test[];
//     for (int i = 0; i <= 1; i++)
//       if (!is_refined(neighbor(2*i-1)) && neighbor(2*i-1).neighbors &&
// 	  (is_local(cell) || is_local(neighbor(2*i-1))))
// 	for (int j = 0; j <= 1; j++)
// 	  for (int k = 0; k <= 1; k++)
// 	    fine(v.x,2*i,j,k) = 1;//v.x[i];
// }
// static inline void face_average (Point point, vector v)
// {
//   foreach_dimension() {
//     #if dimension == 1
//       v.x[] = fine(v.x,0);
//       v.x[1] = fine(v.x,2);
//     #elif dimension == 2
//       v.x[] = (fine(v.x,0,0) + fine(v.x,0,1))/2.;
//       v.x[1] = (fine(v.x,2,0) + fine(v.x,2,1))/2.;
//     #else // dimension == 3
//       v.x[] = (fine(v.x,0,0,0) + fine(v.x,0,1,0) +
// 	       fine(v.x,0,0,1) + fine(v.x,0,1,1))/4.;
//       v.x[1] = (fine(v.x,2,0,0) + fine(v.x,2,1,0) +
// 		fine(v.x,2,0,1) + fine(v.x,2,1,1))/4.;
//     #endif
//   }
// }
  
// static inline void restriction_face (Point point, scalar s)
// {
//   face_average (point, s.v);
// }


void refine_embed_linear_css_test3 (Point point, scalar s)
//void refine_embed_linear2 (Point point, scalar s) //1
{
  foreach_child() {
    if (!css_test3[]){
       s[] = Tsat00;
    }else {
      assert (coarse(css_test3));
      int i = (child.x + 1)/2, j = (child.y + 1)/2;
#if dimension == 2
      if (coarse(fss_test3.x,i) && coarse(fss_test3.y,0,j) &&
	  (coarse(css_test3) == 1. || coarse(css_test3,child.x) == 1. ||
	   coarse(css_test3,0,child.y) == 1. || coarse(css_test3,child.x,child.y) == 1.)) {
	assert (coarse(css_test3,child.x) && coarse(css_test3,0,child.y));
	if (coarse(fss_test3.x,i,child.y) && coarse(fss_test3.y,child.x,j)) {
	  // bilinear interpolation
	  assert (coarse(css_test3,child.x,child.y));
	  s[] = (9.*coarse(s) + 
		 3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
		 coarse(s,child.x,child.y))/16.;
	}
	else
	  // triangular interpolation	  
	  s[] = (2.*coarse(s) + coarse(s,child.x) + coarse(s,0,child.y))/4.;
      }
      else if (coarse(css_test3,child.x,child.y) &&
	       ((coarse(fss_test3.x,i) && coarse(fss_test3.y,child.x,j)) ||
		(coarse(fss_test3.y,0,j) && coarse(fss_test3.x,i,child.y)))) {
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y))/4.;
      }
#else // dimension == 3
      int k = (child.z + 1)/2;
      if (coarse(fss_test3.x,i) > 0.25 && coarse(fss_test3.y,0,j) > 0.25 &&
	  coarse(fss_test3.z,0,0,k) > 0.25 &&
	  (coarse(css_test3) == 1. || coarse(css_test3,child.x) == 1. ||
	   coarse(css_test3,0,child.y) == 1. || coarse(css_test3,child.x,child.y) == 1. ||
	   coarse(css_test3,0,0,child.z) == 1. || coarse(css_test3,child.x,0,child.z) == 1. ||
	   coarse(css_test3,0,child.y,child.z) == 1. ||
	   coarse(css_test3,child.x,child.y,child.z) == 1.)) {
	assert (coarse(css_test3,child.x) && coarse(css_test3,0,child.y) &&
		coarse(css_test3,0,0,child.z));
	if (coarse(fss_test3.x,i,child.y) && coarse(fss_test3.y,child.x,j) &&
	    coarse(fss_test3.z,child.x,child.y,k) &&
	    coarse(fss_test3.z,child.x,0,k) && coarse(fss_test3.z,0,child.y,k)) {
	  assert (coarse(css_test3,child.x,child.y) && coarse(css_test3,child.x,0,child.z) &&
		  coarse(css_test3,0,child.y,child.z) &&
		  coarse(css_test3,child.x,child.y,child.z));
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
      else if (coarse(css_test3,child.x,child.y,child.z) &&
	       ((coarse(fss_test3.z,child.x,child.y,k) &&
		 ((coarse(fss_test3.x,i) && coarse(fss_test3.y,child.x,j)) ||
		  (coarse(fss_test3.y,0,j) && coarse(fss_test3.x,i,child.y))))
		||
		(coarse(fss_test3.z,0,0,k) &&
		 ((coarse(fss_test3.x,i,0,child.z) && coarse(fss_test3.y,child.x,j,child.z)) ||
		  (coarse(fss_test3.y,0,j,child.z) && coarse(fss_test3.x,i,child.y,child.z))))
		||
		(coarse(fss_test3.z,child.x,0,k) &&
		 coarse(fss_test3.x,i) && coarse(fss_test3.y,child.x,j,child.z))
		||
		(coarse(fss_test3.z,0,child.y,k) &&
		 coarse(fss_test3.y,0,j) && coarse(fss_test3.x,i,child.y,child.z))
		))
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y,child.z))/4.;
#endif // dimension == 3
      else {
	// Pathological cases, use 1D gradients.
	s[] = coarse(s);
	foreach_dimension() {
	  if (coarse(fss_test3.x,(child.x + 1)/2) && coarse(css_test3,child.x))
	    s[] += (coarse(s,child.x) - coarse(s))/4.;
	  else if (coarse(fss_test3.x,(- child.x + 1)/2) && coarse(css_test3,- child.x))
	    s[] -= (coarse(s,- child.x) - coarse(s))/4.;
	}
      }
    }
  } 
}

void refine_embed_linear_css_test3_n (Point point, scalar s)
//void refine_embed_linear2 (Point point, scalar s) //1
{
  foreach_child() {
    if (!css_test3_n[]){
     //s[] = 0.;
      if(s.boundary[0] != s.boundary_homogeneous[0])
       s[] = Tsat00;
     else
       s[] = 0.0;
    }else {
      assert (coarse(css_test3_n));
      int i = (child.x + 1)/2, j = (child.y + 1)/2;
#if dimension == 2
      if (coarse(fss_test3_n.x,i) && coarse(fss_test3_n.y,0,j) &&
	  (coarse(css_test3_n) == 1. || coarse(css_test3_n,child.x) == 1. ||
	   coarse(css_test3_n,0,child.y) == 1. || coarse(css_test3_n,child.x,child.y) == 1.)) {
	assert (coarse(css_test3_n,child.x) && coarse(css_test3_n,0,child.y));
	if (coarse(fss_test3_n.x,i,child.y) && coarse(fss_test3_n.y,child.x,j)) {
	  // bilinear interpolation
	  assert (coarse(css_test3_n,child.x,child.y));
	  s[] = (9.*coarse(s) + 
		 3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
		 coarse(s,child.x,child.y))/16.;
	}
	else
	  // triangular interpolation	  
	  s[] = (2.*coarse(s) + coarse(s,child.x) + coarse(s,0,child.y))/4.;
      }
      else if (coarse(css_test3_n,child.x,child.y) &&
	       ((coarse(fss_test3_n.x,i) && coarse(fss_test3_n.y,child.x,j)) ||
		(coarse(fss_test3_n.y,0,j) && coarse(fss_test3_n.x,i,child.y)))) {
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y))/4.;
      }
#else // dimension == 3
      int k = (child.z + 1)/2;
      if (coarse(fss_test3_n.x,i) > 0.25 && coarse(fss_test3_n.y,0,j) > 0.25 &&
	  coarse(fss_test3_n.z,0,0,k) > 0.25 &&
	  (coarse(css_test3_n) == 1. || coarse(css_test3_n,child.x) == 1. ||
	   coarse(css_test3_n,0,child.y) == 1. || coarse(css_test3_n,child.x,child.y) == 1. ||
	   coarse(css_test3_n,0,0,child.z) == 1. || coarse(css_test3_n,child.x,0,child.z) == 1. ||
	   coarse(css_test3_n,0,child.y,child.z) == 1. ||
	   coarse(css_test3_n,child.x,child.y,child.z) == 1.)) {
	assert (coarse(css_test3_n,child.x) && coarse(css_test3_n,0,child.y) &&
		coarse(css_test3_n,0,0,child.z));
	if (coarse(fss_test3_n.x,i,child.y) && coarse(fss_test3_n.y,child.x,j) &&
	    coarse(fss_test3_n.z,child.x,child.y,k) &&
	    coarse(fss_test3_n.z,child.x,0,k) && coarse(fss_test3_n.z,0,child.y,k)) {
	  assert (coarse(css_test3_n,child.x,child.y) && coarse(css_test3_n,child.x,0,child.z) &&
		  coarse(css_test3_n,0,child.y,child.z) &&
		  coarse(css_test3_n,child.x,child.y,child.z));
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
      else if (coarse(css_test3_n,child.x,child.y,child.z) &&
	       ((coarse(fss_test3_n.z,child.x,child.y,k) &&
		 ((coarse(fss_test3_n.x,i) && coarse(fss_test3_n.y,child.x,j)) ||
		  (coarse(fss_test3_n.y,0,j) && coarse(fss_test3_n.x,i,child.y))))
		||
		(coarse(fss_test3_n.z,0,0,k) &&
		 ((coarse(fss_test3_n.x,i,0,child.z) && coarse(fss_test3_n.y,child.x,j,child.z)) ||
		  (coarse(fss_test3_n.y,0,j,child.z) && coarse(fss_test3_n.x,i,child.y,child.z))))
		||
		(coarse(fss_test3_n.z,child.x,0,k) &&
		 coarse(fss_test3_n.x,i) && coarse(fss_test3_n.y,child.x,j,child.z))
		||
		(coarse(fss_test3_n.z,0,child.y,k) &&
		 coarse(fss_test3_n.y,0,j) && coarse(fss_test3_n.x,i,child.y,child.z))
		))
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y,child.z))/4.;
#endif // dimension == 3
      else {
	// Pathological cases, use 1D gradients.
	s[] = coarse(s);
	foreach_dimension() {
	  if (coarse(fss_test3_n.x,(child.x + 1)/2) && coarse(css_test3_n,child.x))
	    s[] += (coarse(s,child.x) - coarse(s))/4.;
	  else if (coarse(fss_test3_n.x,(- child.x + 1)/2) && coarse(css_test3_n,- child.x))
	    s[] -= (coarse(s,- child.x) - coarse(s))/4.;
	}
      }
    }
  } 
}

void restriction_embed_linear_css_test3 (Point point, scalar s)
//void restriction_embed_linear2 (Point point, scalar s)
{  
  // 0 children
  if (!css_test3[] ) {
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
      if (fine(css_test3,0,i,j) && fine(css_test3,1,!i,!j))
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
    if (css_test3[])
      p.x += x, p.y += y, p.z += z, val += s[], nv++;
  assert (nv > 0.);
  s[] = val/nv;
}