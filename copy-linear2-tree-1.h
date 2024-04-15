#include "heights.h"
#include "fractions.h"
extern scalar ff;
extern scalar T;
extern scalar Tlff,Tgff;
extern scalar fs_solid;
extern double Tsat0;
extern double Tsat00;
extern scalar T_solid;
#define EPS 0.0000000001
#define EPS4 0.0000000001
#define EPS3 0.0000000001
#define is_phase(v) (v>=0.5)

#include "./graddiriT.h" 
//#include "./vof2front-advanced.h"
#include "./vof2front-advanced-copy.h"

extern vector hhh;
extern face vector modphase0,modphase1;
// // // extern face vector DD;

extern coord mpi_coord,mpi_dimen;

//for mov_interface
extern scalar deltac;
extern int globali,outstep;
extern bool out_flag;

//for gradient
extern int level_interface;

//#if _MPI
//@def is_boundary_box2(point) (((point.i < GHOSTS) && (mpi_coord.x==0)) || 
//                              ((point.i >= point.n + GHOSTS) && (mpi_coord.x==(mpi_dimen.x-1))) ||
//			                        ((point.j < GHOSTS) && (mpi_coord.y==0)) || 
//                              ((point.j >= point.n + GHOSTS) && (mpi_coord.y==(mpi_dimen.y-1))))
//@
//#else
 #define  is_boundary_box2(cell) is_boundary(cell)

//#endif


#if TREE

attribute {
   scalar ff2;
   double aa;
   double bb;
}
//for temperature restirction
// similar to fraction.h: fraction_refine
void energy_refine(Point point, scalar s){
    scalar ff2 = s.ff2;
    double aa = s.aa;
    double bb = s.bb;
    //scalar n;
    double cc = ff2[];
    double ss = s[];
    if(cc <= 0.0 || cc>=1.0){
        foreach_child(){
           s[] = ss;
        }
    }else{
        coord n = mycs(point,ff2);
        double alpha = plane_alpha (cc, n);
        foreach_child(){
            static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
            coord nc;
            foreach_dimension(){
                 nc.x = child.x * n.x;
            }
            double cc2;
            cc2 = rectangle_fraction(nc,alpha,a,b);
            //s[] = cc2*Trhol*Tcpl + (1-cc2)*Trhog*Tcpg;
             s[] = (ss*cc2*aa + ss*(1-cc2)*bb)/(cc2*aa + (1-cc2)*bb);
        }
    }
}

void energy_refine2(Point point, scalar s){
    scalar ff2 = s.ff2;
    double aa = s.aa;
    double bb = s.bb;
    //scalar n;
    double cc = ff2[];
    double ss = s[];
    if(cc <= 0.0 || cc>=1.0){
        foreach_child(){
           s[] = ss;
        }
    }else{
        coord n = mycs(point,ff2);
        double alpha = plane_alpha (cc, n);
        foreach_child(){
            static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
            coord nc;
            foreach_dimension(){
                 nc.x = child.x * n.x;
            }
            double cc2;
            cc2 = rectangle_fraction(nc,alpha,a,b);
            //s[] = cc2*Trhol*Tcpl + (1-cc2)*Trhog*Tcpg;

           if(cc2 >= 0.5){
              double sum=0.0;
              double weight = 0.0;
              int lim_x = 2*child.x;
              int lim_y = 2*child.y;
              for(int ii=0; abs(ii)<abs(2*child.x);ii+=child.x){
                   for(int jj=0; abs(jj)<abs(2*child.y);jj+=child.y){
                        if(coarse(ff2,ii,jj)>=0.5){
                            sum = sum + coarse(ff2,ii,jj)*coarse(s,ii,jj);
                            weight = weight + coarse(ff2,ii,jj);
                        }
                   }
              }
              if(weight>1e-12){
                    s[] = sum/weight;
              }else{
                   if(s.boundary[0] != s.boundary_homogeneous[0])
                       s[] = Tsat00;
                   else
                       s[] = 0.0;
              }
           }else{
             if(s.boundary[0] != s.boundary_homogeneous[0])
               s[] = Tsat00;
             else
               s[] = 0.0;
           }
        }
    }
}

void energy_refine22(Point point, scalar s){
    scalar ff2 = s.ff2;
    double aa = s.aa;
    double bb = s.bb;
    //scalar n;
    double cc = ff2[];
    double ss = s[];
    if(cc <= 0.0 || cc>=1.0){
        foreach_child(){
           s[] = ss;
        }
    }else{
        coord n = mycs(point,ff2);
        double alpha = plane_alpha (cc, n);
        foreach_child(){
            static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
            coord nc;
            foreach_dimension(){
                 nc.x = child.x * n.x;
            }
            double cc2;
            cc2 = rectangle_fraction(nc,alpha,a,b);
            //s[] = cc2*Trhol*Tcpl + (1-cc2)*Trhog*Tcpg;
           coord p_child;
           p_child.x = child.x *0.25;
           p_child.y = child.y *0.25;
           bool too_close1=false;
           if(cc2 >= 0.5){
              double sum=0.0;
              double weight = 0.0;
              int lim_x = 2*child.x;
              int lim_y = 2*child.y;
              for(int ii=0; abs(ii)<abs(2*child.x);ii+=child.x){
                   for(int jj=0; abs(jj)<abs(2*child.y);jj+=child.y){
                        if(coarse(ff2,ii,jj)>=0.5 && (!too_close1)){
                            //method1-volume-average
                            //sum = sum + coarse(ff2,ii,jj)*coarse(s,ii,jj);
                            //weight = weight + coarse(ff2,ii,jj);
                            //method2-distance-based
                            coord big_nei;
                            big_nei.x = ii;
                            big_nei.y = jj;
                            double dd = sqrt(sq(big_nei.x - p_child.x)+sq(big_nei.y-p_child.y));
                            if(dd<1e-14){
                                 s[] = coarse(s,ii,jj);
                                 too_close1 = true;
                            }else{
                                 weight  = weight +  1.0/dd;
                                 sum = sum + (1.0/dd)*coarse(s,ii,jj);
                            }
                        }
                   }
              }
              if(!too_close1){
                  if(weight>1e-14){
                        s[] = sum/weight;
                  }else{
                      if(s.boundary[0] != s.boundary_homogeneous[0])
                              s[] = Tsat00;
                      else
                              s[] = 0.0;
                  }
              }
           }else{
             if(s.boundary[0] != s.boundary_homogeneous[0])
               s[] = Tsat00;
             else
               s[] = 0.0;
           }
        }
    }
}




void restriction_cprho_average (Point point, scalar s)
{
  double sum = 0.;
  foreach_child()
    sum += cm[]*s[];
  s[] = sum/(1 << dimension)/(cm[] + 1e-30);
}

void restriction_cprho_average2 (Point point, scalar s)
{
  double sum = 0.;
  int num = 0;
  double weight=0.0;
   scalar ff2 = s.ff2;
  if(ff2[]<0.5){
      if(s.boundary[0] != s.boundary_homogeneous[0])
               s[] = Tsat00;
      else
               s[] = 0.0;

      return ;
  }
  foreach_child()
      if(ff2[]>=0.5){
         sum += ff2[]*s[];
         num++;
         weight += ff2[];
       }
  if(num>0)
    s[] = sum/(weight + 1e-30);
  else{
    if(s.boundary[0] != s.boundary_homogeneous[0])
               s[] = Tsat00;
    else
               s[] = 0.0;
  }
}

void restriction_cprho_average23 (Point point, scalar s)
{
  double sum = 0.;
  int num = 0;
  double weight=0.0;
   scalar ff2 = s.ff2;
   bool flag_Tsat0_0 = s.flag_Tsat0_0;
  // if(ff2[]<0.5){
  //     if(s.boundary[0] != s.boundary_homogeneous[0])
  //              s[] = Tsat00;
  //     else
  //              s[] = 0.0;

  //     return ;
  // }
  // foreach_child(){
  //     if(ff2[]>=0.5){
  //        sum += ff2[]*s[];
  //        num++;
  //        weight += ff2[];
  //      }
  // }
  foreach_child(){
      //if(ff2[]>=0.5){
        if(ff2[]>EPS){
         sum += ff2[]*s[];
         num++;
         weight += ff2[];
       }
  }
  if(num>0)
    s[] = sum/(weight + 1e-30);
  else{
    // if(s.boundary[0] != s.boundary_homogeneous[0])
    //            s[] = Tsat00;
    // else
    //            s[] = 0.0;
    if(!flag_Tsat0_0)
               s[] = Tsat00;
    else
               s[] = 0.0;
  }
}

void restriction_cprho_average24 (Point point, scalar s)
{
  double sum = 0.;
  int num = 0;
  double weight=0.0;
   scalar ff2 = s.ff2;
   bool flag_Tsat0_0 = s.flag_Tsat0_0;
  // if(ff2[]<0.5){
  //     if(s.boundary[0] != s.boundary_homogeneous[0])
  //              s[] = Tsat00;
  //     else
  //              s[] = 0.0;

  //     return ;
  // }
  // foreach_child(){
  //     if(ff2[]>=0.5){
  //        sum += ff2[]*s[];
  //        num++;
  //        weight += ff2[];
  //      }
  // }
  foreach_child(){
      if(ff2[]>=0.5){
        //if(ff2[]>EPS){
         sum += ff2[]*s[];
         num++;
         weight += ff2[];
       }
  }
  if(num>0)
    s[] = sum/(weight + 1e-30);
  else{
    // if(s.boundary[0] != s.boundary_homogeneous[0])
    //            s[] = Tsat00;
    // else
    //            s[] = 0.0;
    if(!flag_Tsat0_0)
               s[] = Tsat00;
    else
               s[] = 0.0;
  }
}

void Tl_restriction3_T(Point point, scalar s){
//restriction no interface effecct
  scalar ff3 = s.ff3;
  //double cc = ff3[];
  double sum = 0.;
  double cc4;
 // if(!s.inverse){
    cc4 = ff3[];
  //}else{
  //  cc4 = 1.0-ff3[];
 // }
 // int part_cell=0;
    double weight_i=0.0;
    foreach_child(){
            double cc3;
            //if(!s.inverse){
                cc3 = ff3[]; 
            //}else{
            //    cc3 = 1.0 -ff3[];
           // }
	    if(cc3>EPS){
	     //part_cell++;
	    // sum += cc3*s[];  //since
             sum += s[];
             weight_i += cc3 ; 
	    }
    }
     if(weight_i<=1e-6){

            if(s.boundary[0] != s.boundary_homogeneous[0])
               s[] = Tsat00;
            else
               s[] = 0.0;
     }else{
	  s[] = sum/weight_i;
     }

     if(s[]<0){

      printf("restriction<0 s[]=%g\n",s[]);
     }
}


void restriction_cprho_average22 (Point point, scalar s)
{
  double sum = 0.;
  int num = 0;
  double weight=0.0;
   scalar ff2 = s.ff2;
   bool flag1 = false;
  if(ff2[]<0.5){
      if(s.boundary[0] != s.boundary_homogeneous[0]){
               s[] = Tsat00;
      }else{
               s[] = 0.0;
      }
      flag1 = true;
  //    return ;
  }
  if(!flag1){
    bool too_close1 = false;
    double too_close1_value;
    foreach_child(){
      if(!too_close1){
        if(ff2[]>=0.5){
          coord position3;
          position3.x = child.x*0.25;
          position3.y = child.y*0.25;
          double dd = sqrt(sq(position3.x -0)+sq(position3.x-0));
          if(dd<1e-14){
               too_close1 = true;
               too_close1_value = s[];
          }else{
              //sum += ff2[]*s[];
              //num++;
              //weight += ff2[];
              weight += 1.0/dd;
              sum += (weight*s[]);
              num++;
          }
        }
      }
    }
    if(!too_close1){
        if(num>0){
          s[] = sum/(weight + 1e-30);
        }else{
          if(s.boundary[0] != s.boundary_homogeneous[0])
                    s[] = Tsat00;
          else
                    s[] = 0.0;
        }
     }else{
         s[]  = too_close1_value;

     }
  }
}





#include "./linear2-tree-2.h"
#endif

extern scalar css_test;
extern face vector fss_test;

//refine taken from embed.h
static inline void refine_embed_linear2 (Point point, scalar s)
//void refine_embed_linear2 (Point point, scalar s) //1
{
  foreach_child() {
    if (!css_test[]){
     //s[] = 0.;
      if(s.boundary[0] != s.boundary_homogeneous[0])
       s[] = Tsat00;
     else
       s[] = 0.0;
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


static inline
coord embed_gradient2 (Point point, vector u, coord p, double vb ,coord n)
{
  coord dudn;
  foreach_dimension() {
   // bool dirichlet;
   // double vb = u.x.boundary[embed] (point, point, u.x, &dirichlet);
   // if (dirichlet) {
      double val;
      dudn.x = dirichlet_gradient2 (point, u.x, css_test, n, p, vb, &val,fss_test);
   // }
   // else // Neumann
   //   dudn.x = vb;
    if (dudn.x == nodata)
      dudn.x = 0.;
  }
  return dudn;
}

attribute {
  void (* embed_gradient2) (Point, scalar, coord *, double);
}

static inline void restriction_embed_linear_Tsat0 (Point point, scalar s)
//void restriction_embed_linear_Tsat0 (Point point, scalar s) //5
{  
  // 0 children
  if (!css_test[] ) {
     if(s.boundary[0] != s.boundary_homogeneous[0])
       s[] = Tsat00;
     else
       s[] = 0.0;
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

  /**
  If the gradient is defined and if the variable is not using
  homogeneous boundary conditions, we improve the interpolation using
  this information. */
  
  //if (s.embed_gradient && s.boundary[0] != s.boundary_homogeneous[0]) {
  if (s.boundary[0] != s.boundary_homogeneous[0]) {
    coord o = {x,y,z}, g;
    //s.embed_gradient (point, s, &g);
    double vb = Tsat00;
    s.embed_gradient2 (point, s, &g, vb);
    foreach_dimension()
      s[] += (o.x - p.x/nv)*g.x;
  }
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



// foreach_dimension()
// void embed_face_fraction_refine2_x (Point point, scalar s)
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
//////////////////////////////////////////////
//////////////////////////////////////add here
//////////////////////////////////////////////
//for s field restriction based on  css_test
static inline void restriction_embed_linear2 (Point point, scalar s)
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


static inline void restriction_embed_linear2_T (Point point, scalar s)
//void restriction_embed_linear2_T (Point point, scalar s)
{  
  // 0 children
  bool flag_Tsat0_0 = s.flag_Tsat0_0 ;
  // if (!css_test[] ) {
  //      s[] = Tsat00;
  //   return;
  // }
  if(!css_test[]){
      if(flag_Tsat0_0){
          s[] = 0.0;
      }else{
          s[] = Tsat00;
      }
    return ; 
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
//ffor s field restriction based on  css_test2
extern scalar css_test2;
extern face vector fss_test2;
static inline void restriction_embed_linear3 (Point point, scalar s)
//void restriction_embed_linear3 (Point point, scalar s)
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


//refine for filed for css_test2
static inline void refine_embed_linear3 (Point point, scalar s)
//void refine_embed_linear3 (Point point, scalar s)
{
  foreach_child() {
    if (!css_test2[]){
     //s[] = 0.;
      if(s.boundary[0] != s.boundary_homogeneous[0])
       s[] = Tsat00;
     else
       s[] = 0.0;
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


//refine for fss_test2


// foreach_dimension()
// void embed_face_fraction_refine3_x (Point point, scalar s)
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




//extern double Trhol,Trhog,Tcpl,Tcpg;
extern double Tkg,Tkl,Trhog,Trhol,Tcpg,Tcpl,hfg;
extern int Tl_tree_choice;
event init(t=0) {

 //initial heights, as contact.h
 // for (scalar c in interfaces)
 //    if (ff.height.x.i){
 //     heights (ff, ff.height);
 //        fprintf(stderr,"print something\n");
 //    }


  
#if TREE

// foreach_dimension()
//             fss_test.x.refine=fss_test.x.prolongation = embed_face_fraction_refine2_x;
  

// for(int phase_flag4=0;phase_flag4<=1;phase_flag4++){

//       //for embedflux
//          if(phase_flag4 == 1){
// 		 foreach(){
// 		    ff[]=ff_temp2[];
// 		 }     
//          }else{
// 		 foreach(){
// 		    ff[]=1.0-ff_temp2[];
//                   }

//          }
//    boundary({ff});
//    vof2dist(ff, phi20); 
//    fractions (phi20, css_test, fss_test);
//    boundary({css_test,fss_test});


// #if TREE
//   css_test.refine = embed_fraction_refine2;

//   css_test.prolongation = embed_fraction_refine2; //fraction_refine; // this is in fraction
//   foreach_dimension()
//     fss_test.x.prolongation = embed_face_fraction_refine2_x;
// #endif
//    restriction({css_test,fss_test});


// if(phase_flag4 == 1){
//     T.aa = Trhol*Tcpl;
//     T.bb = Trhog*Tcpg;
//     T.ff2 = css_test;
// }else{
//    T.bb = Trhol*Tcpl;
//    T.aa = Trhog*Tcpg;
//    T.ff2 = css_test;
// }

// T.refine = T.prolongation = energy_refine22 ;
// T.restriction = restriction_cprho_average2;

// }





// T.aa = Trhol*Tcpl;
// T.bb = Trhog*Tcpg;
// T.ff2 = ff;

// T.refine = T.prolongation = energy_refine22 ;
// T.restriction = restriction_cprho_average2 ;

//T.refine = T.prolongation = energy_refine ;
//T.restriction = restriction_cprho_average ;

//T.refine = T.prolongation = refine_embed_linear2;
//T.restriction = restriction_embed_linear_Tsat0;


// T.dirty = true;
#endif //TREE


#if TREE
//   if(Tl_tree_choice==1){
// 	T_tree_function(Tlff,ff,Tsat0);
// 	T_tree_function(Tgff,ff,Tsat0);
//   }else if(Tl_tree_choice==0){
//        // T_tree_function2(Tlff,ff);
// 	//T_tree_function2(Tgff,ff);

// /*
//         Tlff.ff2 = ff;
//         Tlff.Trefine = Tlff.prolongation = energy_refine3 ;
//         Tlff.restriction = restriction_cprho_average3 ;

//         Tgff.ff2 = ff;
//         Tgff.Trefine = Tgff.prolongation = energy_refine3 ;
//         Tgff.restriction = restriction_cprho_average3 ;
// */
//   }else if(Tl_tree_choice == 4){
//         vertex scalar phi20[];
//         vof2dist(ff, phi20); 
//         fractions (phi20, css_test, fss_test);
//         scalar ff_temp[];
//         foreach(){
//            ff_temp[] = ff[];
//            ff[] = 1.0 -ff_temp[];
//         }
//         vof2dist(ff, phi20); 
//         fractions (phi20, css_test2, fss_test2);
//         foreach(){
//            ff[] = ff_temp[];
//         }
//         css_test.refine = embed_fraction_refine2;
//         css_test.prolongation = fraction_refine;
//         foreach_dimension()
//             fss_test.x.prolongation = embed_face_fraction_refine2_x;
//         Tlff.restriction = restriction_embed_linear2;
//         Tlff.refine = Tlff.prolongation = refine_embed_linear2;

//         css_test2.refine = embed_fraction_refine3;
//         css_test2.prolongation = fraction_refine;
//         foreach_dimension()
//             fss_test2.x.prolongation = embed_face_fraction_refine3_x;
//         Tgff.restriction = restriction_embed_linear3;
//         Tgff.refine = Tgff.prolongation = refine_embed_linear3;
//   }
    //  T.restriction = restriction_three_phase1;
    //  T.refine = T.prolongation = refine_three_phase1; // this should be consistent with function used in posson3.h; 
#endif
  
}

/*
event output_height(i++){
      
 vector hh[];
   
 //if(ff.height.x.i){
         heights(ff,hh);
  
	  for(int l=1;l<=depth();l++){
	    char name100[80];
	    sprintf(name100,"height_level%d_%g.dat",l, t);
	    FILE * fp100 = fopen(name100,"a");
	    foreach_level(l){ 
		foreach_dimension(1){
		    int orien = orientation(hh.x[]);
		    if(hh.x[] == nodata) {
		        
		    }else{
		       fprintf (fp100, "%g %g %g %d %g %d %g\n", x, y, ff[], l, hh.x[], orien,height(hh.x[]));
		    }
		 }
	    }
	    fclose(fp100);
	  }
   //  }  


}
*/


#define EPS5 0.000000001
#define CCCMIN CCMIN
void get_T_from_Tlff_Tgff3(int phase_flag,scalar ff6){ //from temperature of fraction of liquid and gass to  temperature of center;

// this should be get from each level??????????

///clamp(ff,0,1)  //////////////////////////////////////////////////////////////////////>??????????????????????????????????????
//?????????????????????????????????????????
  //  foreach(){
  //             if(phase_flag==1){ //solving Tl equation
  //                   if(ff[]>EPS){
  //                     T[] = Tlff[]/ff[];
  //                   }else{
  //                     T[] = Tsat00;
  //                   }
  //             }else{  //solving Tg equation
  //                   if((1-ff[])>EPS){
  //                     T[] = Tgff[]/(1.0-ff[]);
  //                   }else{
  //                     T[] = Tsat00;
  //                   }
  //             }
  //  }
     foreach(){
              if(phase_flag==1){ //solving Tl equation
                    //if(ff6[]>EPS5){
                   if(ff6[]>CCCMIN){
                      T[] = Tlff[];
                      // if(T[]<Tsat00){
                      //     //T[] = Tsat00;
                      //     double tot=0.0;
                      //     double weight=0.0;
                      //     foreach_neighbor(1){
                      //       if(ff6[]>CCCMIN && Tlff[]>Tsat00){
                      //           weight += ff6[];
                      //           tot += ff6[]*Tlff[];
                      //       }
                      //     }
                      //     if(weight>CCCMIN){
                      //       T[] = tot/weight;
                      //     }else{
                      //       T[] = Tsat00;
                      //     }
                      // }
                    }else{
                      T[] = Tsat00;
                    }
              }else{  //solving Tg equation
                    //if((1-ff[])>EPS){
                    // if(ff6[]>EPS5){
                   if(ff6[]>CCCMIN){
                      //T[] = Tsat00;
                          T[] = Tgff[]; 
                      //    if(T[]<Tsat00){
                      //     //T[] = Tsat00;
                      //     double tot=0.0;
                      //     double weight=0.0;
                      //     foreach_neighbor(1){
                      //       if(ff6[]>CCCMIN && Tgff[]>Tsat00){
                      //           weight+=ff6[];
                      //           tot += ff6[]*Tgff[];
                      //       }
                      //     }
                      //     if(weight>CCCMIN){
                      //       T[] = tot/weight;
                      //     }else{
                      //       T[] = Tsat00;
                      //     }
                      // }
                    }else{
                      T[] = Tsat00;
                    }
              }
              // if(ff6[]>=0.5 && ff6[]<=0.5+EPS5){
              //         T[] = Tsat00;
              // }
   }
/*
    char name_temp2[150];
    sprintf(name_temp2,"outfacets/non2-mormal-pid%d.dat",pid());
    FILE * fp111 = fopen(name_temp2,"w");
    foreach(){
            if(T[]<Tsat00 && ff6[] >CCCMIN){
                fprintf(fp111,"itime=%d ff=%g Tlff[]=%g its neightbor:\n",globali,ff6[],Tlff[]);
                int ii=0;
                foreach_neighbor(1){
                   ii = ii+1;
                   fprintf(fp111,"   ii=%d,ff=%g,Tlff=%g\n",ii,ff6[],Tlff[]);
                }
            }
        }
    fclose(fp111);
    MPI_Barrier (MPI_COMM_WORLD);
    char command11[150];
    sprintf(command11, "LC_ALL=C cat outfacets/non2-mormal-pid*.dat > outfacets/non2-mormal-%g",t);
          system(command11);

          char command77[150];
          sprintf(command77, "LC_ALL=C rm -rf  outfacets/non2-mormal-pid*.dat");
          system(command77);
*/
}


void  get_Tlff_Tgff_from_T3(int phase_flag,scalar ff6){
      // foreach(){
      //     if(phase_flag==1){
      //               if(ff[]>EPS){
      //                 Tlff[] = T[]*ff[];
      //               }else{
      //                 Tlff[] = Tsat00*ff[];
      //               }
      //     }else{
      //               if((1.0-ff[])>EPS){
      //                 Tgff[] = T[]*(1.0-ff[]);
      //               }else{
      //                 Tgff[] = Tsat00*(1.0-ff[]);
      //               }
      //     }       
      // }
      // // foreach(){
      // //     if(phase_flag==1){
      // //               if(ff[]>EPS){
      // //                 Tlff[] = T[]*ff[];
      // //               }else{
      // //                 Tlff[] = Tsat00*ff[];
      // //               }
      // //     }else{
      // //               if((1.0-ff[])>EPS){
      // //                 Tgff[] = T[]*(1.0-ff[]);
      // //               }else{
      // //                 Tgff[] = Tsat00*(1.0-ff[]);
      // //               }
      // //     }       
      // // }
       foreach(){
          if(phase_flag==1){
                    //if(ff6[]>EPS5){
                    if(ff6[]>CCCMIN){
                      Tlff[] = T[];
                      // if(Tlff[]<Tsat00){
                      //       Tlff[] = Tsat00;
                      // }
                    }else{
                      Tlff[] = Tsat00;
                    }
          }else{
                   // if(ff6[]>EPS5){
                    if(ff6[]>CCCMIN){
                        Tgff[] = T[];
                        // if(Tgff[]<Tsat00){
                        //      Tgff[] = Tsat00;
                        // }
                    }else{
                      Tgff[] = Tsat00;
                    }
          }       
      }
}



void get_T_from_Tlff_Tgff(){ //from temperature of fraction of liquid and gass to  temperature of center;

// this should be get from each level??????????

///clamp(ff,0,1)  //////////////////////////////////////////////////////////////////////>??????????????????????????????????????
//?????????????????????????????????????????
int choice =1; //=1 for Tl_change
     foreach(){

if(choice ==0){
          if(ff[]>1.0 - EPS && ff[]<1.0){
              // Tlff[] = Tsat0 * ff[];
              //T[] = Tsat0;
              T[] = Tlff[]/ff[];
           }else if(ff[]<EPS && ff[]>0.0){
               //Tgff[] = Tsat0*(1.0-ff[]);
               //T[] = Tsat0;
               T[] = Tgff[]/(1-ff[]);
           }else if(ff[]>EPS && ff[]<(1.0-EPS)){
          
              //coord na= interface_normal3 (point, ff,hhh),p1a;
              coord na= interface_normal3 (point, ff,hhh),p1a;
	      double alphaa = plane_alpha (ff[], na);
             
              double lengtha = line_length_center (na, alphaa, &p1a);//2-dimension,p is centriod of the interface
              //p1a.x = p1a.x*Delta, p1a.y= p1a.y*Delta;
              coord pla;
              plane_center(na,alphaa,ff[],&pla); //center of liquid
              //pla.x = pla.x*Delta, pla.y = pla.y*Delta;
              coord pga;
              pga.x = -ff[]*pla.x/(1.0-ff[]);
              pga.y = -ff[]*pla.y/(1.0-ff[]);
              coord nna;
              nna.x = na.x/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              nna.y = na.y/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              double d1a,d2a,d3a;
              d2a = (p1a.x - 0.0)*nna.x +  (p1a.y - 0.0)*nna.y; // directional distance, toward oppsite normal of interface
              d1a = fabs((p1a.x - pla.x)*nna.x + (p1a.y - pla.y)*nna.y); // normal distance from p1 to pl
              d3a = fabs((p1a.x - pga.x)*nna.x + (p1a.y - pga.y)*nna.y); // normal distance from p1 to pg
              if(fabs(d2a)<EPS4){
                  T[] = Tsat0;
              }else if(fabs(d1a)<EPS4 || fabs(d3a)<EPS4){
                   T[] = Tsat0;  // this is wrong; or get from gradient, we do not handle this.
              }else if(!( fabs(nna.x)<EPS4 || fabs(nna.y)<EPS4)){
                //commet 20220812
                  /*
                  if(ff[]>0.5-EPS3 && ff[]<0.5+EPS3){
                        T[] = Tsat0;
                   }else if(ff[]>0.5+EPS3){  //T calculate using Tlff[]
                        T[]= (d2a/d1a)*(Tlff[]/ff[]) + (d1a - d2a)/d1a*Tsat0; 
                   }else if(ff[]<0.5-EPS3){  // T calculate using Tgff[]
                        T[]=((d2a + d3a)/d3a)*Tsat0 - d2a/d3a*(Tgff[]/(1.0-ff[]));
                   }
                   */
                   //if(ff[]>0.5-EPS3){  //T calculate using Tlff[]
                   if(ff[]>=0.5){  //T calculate using Tlff[]
                        T[]= (d2a/d1a)*(Tlff[]/ff[]) + (d1a - d2a)/d1a*Tsat0; 
                   }else{  // T calculate using Tgff[]
                        T[]=((d2a + d3a)/d3a)*Tsat0 - d2a/d3a*(Tgff[]/(1.0-ff[]));
                   }
              }else{
                 // if(ff[]>0.5-EPS3){
                  if(ff[]>=0.5){
                       T[] = Tlff[]/ff[];
                  }else{ 
                       T[] = Tgff[]/(1-ff[]);
                  }
              }
         }else if(ff[]>=1.0){ //pure liquid cell
             T[] = Tlff[]/ff[]; //ff[]=1.0
         }else if(ff[]<=0.0){
             T[] = Tgff[]/(1.0-ff[]); //(1-ff[])=1.0
         }
      }else if(choice ==1 ){
          /*
                if(ff[]>0.5-EPS3){
              T[]= Tlff[]/ff[];
              Tgff[] = Tsat0*ff[];
                }else{
              T[]= Tgff[]/(1-ff[]);
              Tlff[] = Tsat0*(1-ff[]);
                }
          */         
            if(ff[]>0.5){
                     T[]= Tlff[];
                }else{
              T[]= Tgff[];
              //Tlff[] = Tsat0*(1-ff[]);
                }
          
      }//choice
   }
}



void get_T_from_Tlff_Tgff_forTg(int phase_flag, scalar cc){ //from temperature of fraction of liquid and gass to  temperature of center;
//cc should be original ff

// this should be get from each level??????????

///clamp(ff,0,1)  //////////////////////////////////////////////////////////////////////>??????????????????????????????????????
//?????????????????????????????????????????
int choice =0;
     foreach(){
double eps_new = 1e-10; //CCCCMIN;

if(choice ==0){
          if(cc[]>1.0 - eps_new && cc[]<1.0){
              // Tlff[] = Tsat0 * ff[];
              //T[] = Tsat0;
              if(phase_flag==1)
                 T[] = Tlff[]/cc[];
              else
                T[] = Tsat0;
           }else if(cc[]<eps_new && cc[]>0.0){
               //Tgff[] = Tsat0*(1.0-ff[]);
               //T[] = Tsat0;
               if(phase_flag==1)
                  T[] = Tsat0;
               else
                  T[] = Tgff[]/(1-cc[]);
           }else if(cc[]>eps_new && cc[]<(1.0-eps_new)){
          
              //coord na= interface_normal3 (point, ff,hhh),p1a;
              coord na= interface_normal (point, cc),p1a;
	      double alphaa = plane_alpha (cc[], na);
             
              double lengtha = line_length_center (na, alphaa, &p1a);//2-dimension,p is centriod of the interface
              //p1a.x = p1a.x*Delta, p1a.y= p1a.y*Delta;
              coord pla;
              plane_center(na,alphaa,cc[],&pla); //center of liquid
              //pla.x = pla.x*Delta, pla.y = pla.y*Delta;
              coord pga;
              pga.x = -cc[]*pla.x/(1.0-cc[]);
              pga.y = -cc[]*pla.y/(1.0-cc[]);
              coord nna;
              nna.x = na.x/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              nna.y = na.y/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              double d1a,d2a,d3a;
              d2a = (p1a.x - 0.0)*nna.x +  (p1a.y - 0.0)*nna.y; // directional distance, toward oppsite normal of interface
              d1a = fabs((p1a.x - pla.x)*nna.x + (p1a.y - pla.y)*nna.y); // normal distance from p1 to pl
              d3a = fabs((p1a.x - pga.x)*nna.x + (p1a.y - pga.y)*nna.y); // normal distance from p1 to pg
              if(fabs(d2a)<EPS4){
                  T[] = Tsat0;
              }else if(fabs(d1a)<EPS4 || fabs(d3a)<EPS4){   //distance is too small to get a gradient, delete_small is needed
                      T[] = Tsat0;
              }else if(!( fabs(nna.x)<EPS4 || fabs(nna.y)<EPS4)){
                //commet 20220812
                  /*
                  if(ff[]>0.5-EPS3 && ff[]<0.5+EPS3){
                        T[] = Tsat0;
                   }else if(ff[]>0.5+EPS3){  //T calculate using Tlff[]
                        T[]= (d2a/d1a)*(Tlff[]/ff[]) + (d1a - d2a)/d1a*Tsat0; 
                   }else if(ff[]<0.5-EPS3){  // T calculate using Tgff[]
                        T[]=((d2a + d3a)/d3a)*Tsat0 - d2a/d3a*(Tgff[]/(1.0-ff[]));
                   }
                   */
                   //if(ff[]>0.5-EPS3){  //T calculate using Tlff[]
                   //if(ff[]>=0.5){  //T calculate using Tlff[]
                  if(phase_flag == 1){
                        T[]= (d2a/d1a)*(Tlff[]/cc[]) + (d1a - d2a)/d1a*Tsat0; 
                   }else{  // T calculate using Tgff[]
                        T[]=((d2a + d3a)/d3a)*Tsat0 - d2a/d3a*(Tgff[]/(1.0-cc[]));
                   }
              }else{
                 // if(ff[]>0.5-EPS3){
                 // if(ff[]>=0.5){
                  if(phase_flag==1){
                       T[] = Tlff[]/cc[];
                  }else{ 
                       T[] = Tgff[]/(1-cc[]);
                  }
              }
         }else if(cc[]>=1.0){ //pure liquid cell
            if(phase_flag == 1){
                 T[] = Tlff[]/cc[]; //ff[]=1.0
            }else{
                 T[] = Tsat0;
            }
         }else if(cc[]<=0.0){
           if(phase_flag==1){
                  T[] = Tsat0;
           }else{
              T[] = Tgff[]/(1.0-cc[]); //(1-ff[])=1.0
           }
         }
      }else if(choice ==1 ){
          /*
                if(ff[]>0.5-EPS3){
              T[]= Tlff[]/ff[];
              Tgff[] = Tsat0*ff[];
                }else{
              T[]= Tgff[]/(1-ff[]);
              Tlff[] = Tsat0*(1-ff[]);
                }
          */         
              //         if(ff[]>=0.5){
              // T[]= Tlff[]/ff[];
              // //Tgff[] = Tsat0*ff[];
              //   }else{
              // T[]= Tgff[]/(1-ff[]);
              // //Tlff[] = Tsat0*(1-ff[]);
              //   }
          
      }//choice
   }
}


extern int choice_T_adap;
static void T_gradient_dirichlet(scalar phase0Tgrad, scalar phase1Tgrad){
//void T_gradient_dirichlet(scalar phase0Tgrad, scalar phase1Tgrad){ //4
///T_interface_normal_gradient
//scalar phase0Tgrad[],phase1Tgrad[];
scalar ff_temp[];
scalar css_test[];
face vector fss_test[];
vertex scalar phi0[],phi1[]; //?????????????????????

foreach(){
     phase0Tgrad[]=0.0;
     phase1Tgrad[]=0.0;
}

int phase;

int choice_get_T;
if(choice_T_adap == 1){
       choice_get_T = 1;
}else if(choice_T_adap ==2){
       choice_get_T = 4;
}
// condering jun bu junyun
//int choice_get_T = 4; // 1 get_T_from_Tlff_Tgff
                      // 3 get_T_from_Tlff_Tgff3
                      // 4 get_T_from_Tlff_Tgff_forTg

if(choice_get_T==1){
 //get_T_from_Tlff_Tgff(); //based on get_T_from_Tlff_Tgff, modify EPS to CCCCMIN
}

scalar T_temp[];
if(choice_get_T!=4){
  foreach(){
        T_temp[] = T[];
  }
}



for(phase=0;phase<=1;phase++){
//for(phase=1;phase>=0;phase--){
/*
   foreach(){
      int local_phase;
      local_phase = (int)(ff[]>=(0.5-EPS)) ;
       if(local_phase==0){
         T[] = 0.0;
       }

    }
*/
if(choice_get_T==4){
  //get_T_from_Tlff_Tgff_forTg(phase,ff); // ff should be the original ff

#if 0
    foreach(){
        if(phase==1){
            if(ff[]<EPS){
              T[] = Tsat00;
            }else{
              T[] = Tlff[];
            }
        }else{
            if(1.0-ff[]<EPS){
              T[] = Tsat00;
            }else{
              T[] = Tgff[];
            }
        }
    }
#else
    // to consistent with gradient_function in poisson3, here when ff<CCCCMIN, T should set to Tsat00
    // get_tempa_from_Tlff_Tgff_forTg3 in poisson3
    foreach(){
        if(phase==1){
            if(ff[]<CCCMIN){
              T[] = Tsat00;
            }else{
              T[] = Tlff[];
            }
        }else{
            if(1.0-ff[]<CCCMIN){
              T[] = Tsat00;
            }else{
              T[] = Tgff[];
            }
        }
    }
#endif

  
  

}
   if(phase == 1){
         foreach(){
            ff_temp[]=ff[];
         }     
   }else{
         foreach(){
            ff_temp[]=1.0-ff[];
         }

   }
boundary({ff_temp});

foreach(){
  css_test[] = nodata;
}
foreach_face(){
  fss_test.x[] = nodata;
}
   vof2dist(ff_temp, phi0); 
   fractions (phi0, css_test, fss_test);

/*
if(globali%outstep==0){
   if(phase==1){
          char name29[80];
          sprintf(name29,"outfacets/check_phi1%g.dat", t);
          FILE * fp29 = fopen(name29,"a");
          foreach_vertex(){
              fprintf(fp29,"%g %g %g\n",x,y,phi0[]);
          }
          fclose(fp29);

          char name129[80];
          sprintf(name129,"outfacets/check_css1%g.dat", t);
          FILE * fp129 = fopen(name129,"a");
          foreach(){
              fprintf(fp129,"%g %g %g %g\n",x,y,css_test[],ff_temp[]);
          }
          fclose(fp129);

   }else{
        char name29[80];
          sprintf(name29,"outfacets/check_phi0%g.dat", t);
          FILE * fp29 = fopen(name29,"a");
          foreach_vertex(){
              fprintf(fp29,"%g %g %g\n",x,y,phi0[]);
          }
          fclose(fp29);

          char name129[80];
          sprintf(name129,"outfacets/check_css0%g.dat", t);
          FILE * fp129 = fopen(name129,"a");
          foreach(){
              fprintf(fp129,"%g %g %g %g\n",x,y,css_test[],ff_temp[]);
          }
          fclose(fp129);
   }
}
*/
   boundary({css_test,fss_test});

#if TREE
  css_test.refine = embed_fraction_refine2;
//For prolongation we cannot use the same function since the surface fraction 
//field fss_test is not necessarily defined for prolongation cells.
//So we switch back to the default fraction refinement (which is less accurate but only relies on css_test).

  // fixme: could work better with embed_fraction_refine
  // see porous1.tst
  css_test.prolongation = embed_fraction_refine2; //fraction_refine; // this is in fraction
  foreach_dimension()
    fss_test.x.prolongation = embed_face_fraction_refine2_x;
//Note that we do not need to change the refine method since the default refine method calls the prolongation method for each component.

#endif


   restriction({css_test,fss_test});
//restriction({css_test,fss_test});
//restriction({css_test,fss_test});

/*
char name73[80];
sprintf(name73,"ff_temp_phase%d-%g.dat",phase, t);
FILE * fp73 = fopen(name73,"a");
fprintf(fp73,"#i j ff ff_temp\n");
foreach(){
 //double rr=sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5));
 //if(rr<0.01+0.2 && rr>0.2-0.01 )
    if(ff[]<1.0 && ff[]>0)
       if(point.j==BGHOSTS){
         fprintf(fp73,"%d %d %g %g\n",point.i,point.j-1,ff[0,-1],css_test[0,-1]);
         fprintf(fp73,"%d %d %g %g\n",point.i,point.j,ff[],css_test[]);
       }else{
           if(point.j<10)
            fprintf(fp73,"%d %d %g %g\n",point.i,point.j,ff[],css_test[]);
       }
}
fclose(fp73);
*/

/* 
char name48[80];
sprintf(name48,"css_test%d-%g.dat", phase, t);
FILE * fp48 = fopen(name48,"a");
fprintf(fp48,"# I J ff css_test fss_test \n");
foreach(){
    if(phase==0){
      double temp5 = fabs(ff[]+css_test[]-1.0);
      if(temp5>0){
   // fprintf(fp29,"I=%d J=%d ff=%g phase0Tgrad=%g phase1Tgrad=%g \n", \
point.i,point.j,ff[],phase0Tgrad[],phase1Tgrad[]); 
       fprintf(fp48,"%d %d %g %g %g %g %g\n", \
point.i,point.j,x,y,ff[],css_test[],temp5); 
      }
   }else{
       double temp5 =fabs(ff[]-css_test[]);
       if(temp5>0){
   // fprintf(fp29,"I=%d J=%d ff=%g phase0Tgrad=%g phase1Tgrad=%g \n", \
point.i,point.j,ff[],phase0Tgrad[],phase1Tgrad[]); 
       fprintf(fp48,"%d %d %g %g %g %g %g\n ", \
point.i,point.j,x,y,ff[],css_test[],temp5); 

       }
   }
}
fclose(fp48);


char name49[80];
sprintf(name49,"fss_test%d-%g.dat", phase, t);
FILE * fp49 = fopen(name49,"a");
fprintf(fp48,"# x y fss_test \n");
foreach_face(){
   // fprintf(fp29,"I=%d J=%d ff=%g phase0Tgrad=%g phase1Tgrad=%g \n", \
point.i,point.j,ff[],phase0Tgrad[],phase1Tgrad[]); 
     fprintf(fp49,"x=%g y=%g fss_test=%g\n", x,y,fss_test.x[]); 
}
fclose(fp49);
*/
//////////////////////neumann.c
/*
if(phase=1){
    foreach(){
	 css_test[]=ff[];
    }
}else{
    foreach(){
	 css_test[]=1.0-ff[];
    }
}
*/

//char name74[80];
//sprintf(name74,"nn_phase%d-%g.dat",phase, t);
//FILE * fp74 = fopen(name74,"a");
//fprintf(fp74,"#i j ff nnx nny\n");
if(choice_get_T==1){
    foreach(){
        if((is_phase(ff[]) && (phase ==1)) || (!is_phase(ff[]) && (phase ==0))){
            T[] = T_temp[];
        }else{
            T[] = Tsat00;

        }
    }
}else if(choice_get_T==3){
   // get_T_from_Tlff_Tgff3(phase,css_test);
  get_T_from_Tlff_Tgff3(phase,css_test);
}
double limit2 = 1e-6; //EPS; //1e-6;
foreach() {
  
   if(level == level_interface){
      if ((css_test[] > limit2) && (css_test[] < 1.0-limit2)) {

	//scalar s = T; //scalar s = a;
  
        scalar s=T;
	coord n = facet_normal (point, css_test, fss_test), p;
	double alpha = plane_alpha (css_test[], n);
	// double length = line_length_center (n, alpha, &p); //2-dimension,p is the center point of the interface
	// double xx = x + p.x*Delta, yy = y + p.y*Delta;  //?????????this changed the point forever or please check the outfput(fprintf)
	// double nn = sqrt (sq(n.x) + sq(n.y));
	// n.x /= nn, n.y /= nn;
  double length = plane_area_center (n, alpha, &p); //2-dimension,p is the center point of the interface
	double xx = x + p.x*Delta, yy = y + p.y*Delta, zz = z + p.z*Delta;  //?????????this changed the point forever or please check the outfput(fprintf)
	double nn = sqrt (sq(n.x) + sq(n.y) + sq(n.z));
	n.x /= nn, n.y /= nn, n.z /= nn;

     


	//double dsdn = (n.x*(dsdr*cos(theta) - dsdtheta*sin(theta)) +
	//	       n.y*(dsdr*sin(theta) + dsdtheta*cos(theta)));

	//e[] = dsdn - dirichlet_gradient (point, s, cs, n, p, exact (x, y));
        //5 or 6 xingcan??????????????????
        double temp7;
//	printf("Tsat0=%g Tsat===%g/n",Tsat0,Tsat00);
        //double temp4=dirichlet_gradient2 (point, s, css_test, n, p, Tsat0, &temp7,fss_test);
        double temp4=dirichlet_gradient2 (point, s, ff_temp, n, p, Tsat00, &temp7,fss_test);
        //double temp4 = dirichlet_gradient2 (point, s, ff_temp, n, p, Tsat0, &temp7,fss_test);
        if(phase==1){
           phase1Tgrad[]=temp4;      
        }else{
           // temp4 is caculated when gas=1 fluid=0, so it is calculated in the normal of gas, since we take the normal of the liquid as the normal of interface, and dirichlet_gradient2 is caculate in the direction of (fluid=0 gas=1)so gas normal is opposite to interface normal. to get the gradient in the direction of interface normal, we should add an "-"  
           phase0Tgrad[]=-temp4;
        }
     // m=gasT_gradient-liquidT_gradient
/*
#if 1  //5 %g but 4 variables ????????????????????
       fprintf (stderr, "g %g %g %g %g\n",
		x, y, dsdn,
		dirichlet_gradient (point, s, cs, n, p, exact (x, y)));
#endif
*/ 
/*
 if(point.j==BGHOSTS && phase==1){
        // fprintf(fp73,"%d %d %g %g %g %g\n",point.i,point.j-2,ff[0,-2],css_test[0,-2],n.x,n.y);
       //  fprintf(fp73,"%d %d %g %g %g %g\n",point.i,point.j-1,ff[0,-1],css_test[0,-1],n.x,n.y);
         
       //  fprintf(fp73,"%d %d %g %g %g %g\n",point.i,point.j,ff[],css_test[],n.x,n.y);
   //    fprintf(fp74,"%d %d %g %g %g %g %g %g\n",point.i,point.j,ff[],T[]-Tsat0,css_test[],n.x,n.y,phase1Tgrad[]);

	       if(fss_test.x[0,-1]){
		   fprintf(stderr,"ghost fss_test exit\n");
	       }else{
                   fprintf(stderr,"ghost fss_test doesn't exit\n");
               }
       }else{
           // if( point.j>250)
    //       if(point.j<10 && phase==1)
   //         fprintf(fp74,"%d %d %g %g %g %g %g %g\n",point.i,point.j,ff[],T[]-Tsat0,css_test[],n.x,n.y,phase1Tgrad[]);
       }
        
*/



      }
   } // if level_interface
//      else
//	e[] = nodata;

    }


// // // if(globali%outstep==0){
// // //   if(phase==1){
// // //       char name29[80];
// // //       sprintf(name29,"outfacets/check_Tlff%g.dat", t);
// // //       FILE * fp29 = fopen(name29,"a");
// // //       foreach(){
// // //          if(level == level_interface && (fabs(z-L0/2.0)<=Delta/2.0))
// // //             fprintf(fp29,"%g %g %g %d %g %g %g\n",x,y,ff[],phase,T[],phase1Tgrad[]);
// // //       }
// // //       fclose(fp29);
// // //       //exit(1);
// // //   }
// // //   if(phase==0){
// // //       char name39[80];
// // //       sprintf(name39,"outfacets/check_Tgff%g.dat", t);
// // //       FILE * fp39 = fopen(name39,"a");
// // //       foreach(){
// // //         if(level == level_interface && (fabs(z-L0/2.0)<=Delta/2.0))
// // //              fprintf(fp39,"%g %g %g %d %g %g %g\n",x,y,ff[],phase,T[],phase0Tgrad[]);
// // //       }
// // //       fclose(fp39);
// // //       //exit(1);
// // //   }
// // // }


//fclose(fp74);
//exit(1);

  // //  foreach(){
  // //     if((css_test[] > 0) && (css_test[] < 0)){
  // //        if (!(css_test[] > limit2) && (css_test[] < 1.0-limit2)){
  // //           double temp8;
  // //           double val_tot=0.0;
  // //           double weight_i=0.0;
  // //           double weight_tot=0.0;
  // //           int weight_n=0;
  // //           foreach_neighbor(2){
  // //               if(!is_boundary_box2(cell)){
  // //                 if ((css_test[] > limit2) && (css_test[] < 1.0-limit2)){
  // //                   weight_n +=1;
  // //                   weight_i = 1.0;
  // //                   weight_tot += weight_i;
  // //                   if(phase==1){
  // //                       val_tot+=phase1Tgrad[];      
  // //                     }else{
  // //                       val_tot+=phase0Tgrad[];
  // //                     }
                    
  // //                 }
  // //               }
  // //           }
  // //           if(weight_tot>0){
  // //                if(phase==1){
  // //                       phase1Tgrad[]=val_tot/weight_tot;      
  // //                }else{
  // //                       phase0Tgrad[]=val_tot/weight_tot;
  // //                }
  // //           }

  // //        }
  // //     }

  // //  }


}

if(choice_get_T==1){
  foreach(){
        T[] = T_temp[];
  }
}

/*
char name29[80];
sprintf(name29,"gradTwhole%g.dat", t);
FILE * fp29 = fopen(name29,"a");
fprintf(fp29,"# I J ff Phase0Tgrad Phase1grad \n");
//char name53[80];
//sprintf(name53,"gradTwhole2-%g.dat", t);
//FILE * fp53 = fopen(name53,"a");
//fprintf(fp53,"# I J ff Phase0Tgrad Phase1grad \n");
foreach(){
  if(fabs(phase1Tgrad[])>0.0 || fabs(phase0Tgrad[])>0.0){
   // fprintf(fp29,"I=%d J=%d ff=%g phase0Tgrad=%g phase1Tgrad=%g \n", \
point.i,point.j,ff[],phase0Tgrad[],phase1Tgrad[]);
     double temp6;
     temp6 = phase0Tgrad[] - phase1Tgrad[]; 
     fprintf(fp29,"%d %d %g %g %g %g\n", \
point.i,point.j,ff[],phase0Tgrad[],phase1Tgrad[],temp6);
//     if(phase1Tgrad[]>-4.86){
//         fprintf(fp53,"%d %d %g %g %g %g %g %g %g %g %g\n", \
point.i,point.j,x,y,ff[],phase0Tgrad[],phase1Tgrad[],mod_phase1.x[],mod_phase1.x[1]);
//      }
//     if(fabs(x-0.306640625)<1e-5 && fabs(y-0.52832)<1e-5){
 //         fprintf(fp53,"%d %d %g\n", \
point.i,point.j,mod_phase1.x[]);
//     }
  }
}
fclose(fp29);
//fclose(fp53);

*/

//char name30[80];
//sprintf(name30,"gradTinterface%g.dat", t);
//FILE * fp30 = fopen(name30,"a");
//fprintf(fp30,"# I J r theta ff Phase0Tgrad Phase1grad theory0 theory1 error0 error1\n");
/*
foreach(){
   if( ff[]<1.0-EPS && ff[]>EPS ) {

        // coord n = facet_normal (point, css, fss_test), p;
        coord n= interface_normal (point, ff),p;
	double alpha = plane_alpha (ff[], n);
	double length = line_length_center (n, alpha, &p); //2-dimension,p is the center point of the interface
        double x2,y2;
	x2 = x + p.x*Delta - 0.5, y2 = y + p.y*Delta - 0.5;  //?????????this changed the point forever or please check the outfput(fprintf)
	double theta2 = atan2(y2,x2)*180.0/3.1415926, rrr = sqrt(x2*x2 + y2*y2);
      //parabola erro  T=(r/0.2)**2  dT/dr=r/0.04, interface of water and interface is opposite to the r direction, then theory3=-10
          // double theory3 = -10.0; 
      //linear  T=r/0.2; theory3=-dT/dr=-0.5 directional drective of normal direction of interface
          // double theory3 = -5.0;
      //third-order T=(r/0.2)**3; theory3=-dT/dr=3*(r**2)/(0.2**3)
             double theory3 = -15; 
        double error0, error1;
        error0 =  fabs((phase0Tgrad[]-theory3)/theory3);
        error1 =  fabs((phase1Tgrad[]-theory3)/theory3);
   // fprintf(fp29,"I=%d J=%d ff=%g phase0Tgrad=%g phase1Tgrad=%g \n", \
point.i,point.j,ff[],phase0Tgrad[],phase1Tgrad[]); 
     // double radius1= sqrt(pow(x-0.5,2.0)+pow(y-0.5,2.0));
     // double theta2 = atan2(yc, xc)
     //x=x+6.0;
  //   fprintf(fp30,"%d %d %g %g %g %g %g %g %g %g\n", \
point.i,point.j,rrr,theta2,ff[],phase0Tgrad[],phase1Tgrad[],theory3,error0,error1); 

     //fprintf(fp30,"%d %d %g %g %g %g %g %g\n", \
point.i,point.j,rrr,theta2,ff[],phase1Tgrad[],theory3,error1); 
     
   }

}
*/

/*
foreach(){
   if( ff[]<1.0 && ff[]>0.0 ) {
     fprintf(fp30,"x=%g\n",x); 
}
}
*/
//fclose(fp30);
///////////////////////////////



}


void get_T_from_Tlff_Tgff2(){ //from temperature of fraction of liquid and gass to  temperature of center;

// this should be get from each level??????????

///clamp(ff,0,1)  //////////////////////////////////////////////////////////////////////>??????????????????????????????????????
//?????????????????????????????????????????
int choice =1;
     foreach(){

if(choice ==0){
          if(ff[]>1.0 - EPS && ff[]<1.0){
              // Tlff[] = Tsat0 * ff[];
              //T[] = Tsat0;
              T[] = Tlff[]/ff[];
           }else if(ff[]<EPS && ff[]>0.0){
               //Tgff[] = Tsat0*(1.0-ff[]);
               //T[] = Tsat0;
               T[] = Tgff[]/(1-ff[]);
           }else if(ff[]>EPS && ff[]<(1.0-EPS)){
          
              coord na= interface_normal3 (point, ff,hhh),p1a;
	      double alphaa = plane_alpha (ff[], na);
             
              double lengtha = line_length_center (na, alphaa, &p1a);//2-dimension,p is centriod of the interface
              //p1a.x = p1a.x*Delta, p1a.y= p1a.y*Delta;
              coord pla;
              plane_center(na,alphaa,ff[],&pla); //center of liquid
              //pla.x = pla.x*Delta, pla.y = pla.y*Delta;
              coord pga;
              pga.x = -ff[]*pla.x/(1.0-ff[]);
              pga.y = -ff[]*pla.y/(1.0-ff[]);
              coord nna;
              nna.x = na.x/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              nna.y = na.y/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              double d1a,d2a,d3a;
              d2a = (p1a.x - 0.0)*nna.x +  (p1a.y - 0.0)*nna.y; // directional distance, toward oppsite normal of interface
              d1a = fabs((p1a.x - pla.x)*nna.x + (p1a.y - pla.y)*nna.y); // normal distance from p1 to pl
              d3a = fabs((p1a.x - pga.x)*nna.x + (p1a.y - pga.y)*nna.y); // normal distance from p1 to pg
              if(fabs(d2a)<EPS4 || fabs(d1a)<EPS4 || fabs(d3a)<EPS4){
                  T[] = Tsat0;
              }else if(!( fabs(nna.x)<EPS4 || fabs(nna.y)<EPS4)){
                //commet 20220812
                  /*
                  if(ff[]>0.5-EPS3 && ff[]<0.5+EPS3){
                        T[] = Tsat0;
                   }else if(ff[]>0.5+EPS3){  //T calculate using Tlff[]
                        T[]= (d2a/d1a)*(Tlff[]/ff[]) + (d1a - d2a)/d1a*Tsat0; 
                   }else if(ff[]<0.5-EPS3){  // T calculate using Tgff[]
                        T[]=((d2a + d3a)/d3a)*Tsat0 - d2a/d3a*(Tgff[]/(1.0-ff[]));
                   }
                   */
                   if(ff[]>=0.5){  //T calculate using Tlff[]
                        T[]= (d2a/d1a)*(Tlff[]/ff[]) + (d1a - d2a)/d1a*Tsat0; 
                   }else{  // T calculate using Tgff[]
                        T[]=((d2a + d3a)/d3a)*Tsat0 - d2a/d3a*(Tgff[]/(1.0-ff[]));
                   }
              }else{
                  if(ff[]>=0.5){
                       T[] = Tlff[]/ff[];
                  }else{ 
                       T[] = Tgff[]/(1-ff[]);
                  }
              }
         }else if(ff[]>=1.0){ //pure liquid cell
             T[] = Tlff[]/ff[]; //ff[]=1.0
         }else if(ff[]<=0.0){
             T[] = Tgff[]/(1.0-ff[]); //(1-ff[])=1.0
         }
      }else if(choice ==1 ){
 /*
	     if(ff[]>0.5-EPS3){
		T[]= Tlff[]/ff[];
		Tgff[] = Tsat0*ff[];
	     }else{
		T[]= Tgff[]/(1-ff[]);
		Tlff[] = Tsat0*(1-ff[]);
	    }
*/         
            if(ff[]>=0.5){
		T[]= Tlff[]/ff[];
		//Tgff[] = Tsat0*ff[];
	     }else{
		T[]= Tgff[]/(1.0-ff[]);
		//Tlff[] = Tsat0*(1-ff[]);
	    }

            if(ff[]>=0.5 && ff[]<=0.5+EPS){
                T[] = Tsat0;
            }else if(ff[]<=0.5 && ff[]>0.5-EPS){
                T[] = Tsat0;
            }
            
 
      }//choice
   }
}






scalar phase0Tgrad[],phase1Tgrad[];


void  get_Tlff_Tgff_from_T(){  // from T of cell center to bycenter of fluid.  This is from Tsat0 and \partial_T/partian_n
                       // another method is from T_sat T_liquid(for liquid cell) to approximate, for gas cell still using the Tsat0 and partial_T/partian_n

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// in diffusion step, Temperature in mixed cell is temperature in centerpoint temperature.
// in vof step, Temperature becomes temperature in bycenter(zhongxin) temperature. 
// This difference may be an improvement for the temperature diffusion equation.
// in ff>0.5 Tlff = temperature of the cell center point, while Tgff = Temprature of the bycenter(zhongxin)


//scalar phase0Tgrad[],phase1Tgrad[];  // 2022-0907-need to be deleted


    T_gradient_dirichlet(phase0Tgrad,phase1Tgrad);
//char name35[80];
//sprintf(name35,"Tlff_Tgff_mixed%g.dat", t);
//FILE * fp35 = fopen(name35,"a");

//char name56[80];
//sprintf(name56,"gas_cell_liquid_temp%g.dat", t);
//FILE * fp56 = fopen(name56,"a");
//fprintf(fp56,"#i j x y Tlff Tgff\n");
//char name57[80];
//sprintf(name57,"liquid_cell_gas_temp%g.dat", t);
//FILE * fp57 = fopen(name57,"a");
//fprintf(fp57,"#i j x y Tlff Tgff\n");


int choice=1;
  foreach(){
      
  if(choice == 1)   {
      if(ff[]< 0.5-EPS){ //gas cell
           if(ff[]>EPS){  //from center temperature to fraction temperature , 
             ////for temperature of fration greater than 0.5, using T[], Tsat0
              coord na= interface_normal (point, ff),p1a;
	      double alphaa = plane_alpha (ff[], na);
              ////p1a.x = p1a.x*Delta, p1a.y= p1a.y*Delta;
              double lengtha = line_length_center (na, alphaa, &p1a);//2-dimension,p is centriod of the interface
              
              coord pla;
              plane_center(na,alphaa,ff[],&pla); //center of liquid
              ////pla.x = pla.x*Delta, p1a.y= pla.y*Delta;
              coord pga;
              pga.x = -ff[]*pla.x/(1.0-ff[]);
              pga.y = -ff[]*pla.y/(1.0-ff[]);
              coord nna;
              nna.x = na.x/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              nna.y = na.y/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              double d1a,d2a,d3a;
              d2a = (p1a.x - 0.0)*nna.x +  (p1a.y - 0.0)*nna.y; // directional distance, toward oppsite normal of interface
              d1a = fabs((p1a.x - pla.x)*nna.x + (p1a.y - pla.y)*nna.y); // normal distance from p1 to pl
              d3a = fabs((p1a.x - pga.x)*nna.x + (p1a.y - pga.y)*nna.y); // normal distance from p1 to pg
              if(fabs(d2a)<EPS4  || fabs(d1a)<EPS4 || fabs(d3a)<EPS4 || fabs(nna.x)<EPS4 || fabs(nna.y)<EPS4){ // using \partial_T/\partial_n,   T[] from diffusion is not that accurate.
                  coord n= interface_normal3 (point, ff,hhh),p1,p2,p3;
	          double alpha = plane_alpha (ff[], n);
	          double length = line_length_center (n, alpha, &p1);//2-dimension,p is the center point of the interface
                  double x1,y1;
	          x1 = x + p1.x*Delta, y1 = y + p1.y*Delta;  
		  double x2,y2;
		  plane_center(n,alpha,ff[],&p3);  // center of liquid
		  p2.x = -ff[]*p3.x/(1.0-ff[]);
		  p2.y = -ff[]*p3.y/(1.0-ff[]);
		  x2 = x + p2.x*Delta, y2 = y + p2.y*Delta;
		  coord dp1top2,nn;
		  dp1top2.x = x2 - x1;
		  dp1top2.y = y2 - y1;
		  nn.x = n.x/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
		  nn.y = n.y/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
		  double dd= dp1top2.x * nn.x + dp1top2.y * nn.y;
                              
                  Tgff[] = Tsat0 + phase0Tgrad[]*dd;
              }else{
                  Tgff[] = Tsat0*((d2a + d3a)/(d2a)) - T[]*(d3a/(d2a));
              }
            }else if(ff[]>0.0){ 
                Tgff[] = T[];//////// COMMENT 20220812 Tsat0; this comment is error
               // Tgff[] = Tsat0;
            }else{//pure gas cell(ff<)
                Tgff[] = T[];
            }
 
           if(ff[]>EPS){  
              coord n= interface_normal3 (point, ff, hhh),p1,p2;
	      double alpha = plane_alpha (ff[], n);
	      double length = line_length_center (n, alpha, &p1);//2-dimension,p is the center point of the interface
              double x1,y1;
	      x1 = x + p1.x*Delta, y1 = y + p1.y*Delta;  
              double x2,y2;
              plane_center(n,alpha,ff[],&p2);  // center of liquid
              x2 = x + p2.x*Delta, y2 = y + p2.y*Delta;
              coord dp1top2,nn;
              dp1top2.x = x2 - x1;
              dp1top2.y = y2 - y1;
              nn.x = n.x/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
              nn.y = n.y/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
              double dd= dp1top2.x * nn.x + dp1top2.y * nn.y;
                              
              Tlff[] = Tsat0 + phase1Tgrad[]*dd;
              //fprintf(fp56,"%d %d %g %g %g %g\n",point.i,point.j,x,y,Tlff[],Tgff[]);
             // fprintf(fp35,"l: i=%d j=%d ff=%g Tlff=%g Tgff=%g dd=%g nnx=%g nny=%g d12x=%g d12y=%g phase1Tgrad=%g \n",point.i,point.j, ff[], Tlff[], Tgff[], dd, nn.x, nn.y,dp1top2.x, dp1top2.y,phase1Tgrad[]);
 
             }else{
                 Tlff[] = Tsat0;
             }  
             
      }else{  //liquid cell 
           if(ff[]<1.0-EPS){  //from center temperature to fraction temperature , 
            
              ////for temperature of fration greater than 0.5, using T[], Tsat0
              coord na= interface_normal (point, ff),p1a;
	      double alphaa = plane_alpha (ff[], na);
             
              double lengtha = line_length_center (na, alphaa, &p1a);//2-dimension,p is centriod of the interface
              ////p1a.x = p1a.x*Delta, p1a.y= p1a.y*Delta;
              coord pla;
              plane_center(na,alphaa,ff[],&pla); //center of liquid
              ////pla.x = pla.x*Delta, pla.y= pla.y*Delta;
              coord pga;
              pga.x = -ff[]*pla.x/(1.0-ff[]);
              pga.y = -ff[]*pla.y/(1.0-ff[]);
              coord nna;
              nna.x = na.x/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              nna.y = na.y/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              double d1a,d2a,d3a;
              d2a = (p1a.x - 0.0)*nna.x +  (p1a.y - 0.0)*nna.y; // directional distance, toward oppsite normal of interface
              d1a = fabs((p1a.x - pla.x)*nna.x + (p1a.y - pla.y)*nna.y); // normal distance from p1 to pl
              d3a = fabs((p1a.x - pga.x)*nna.x + (p1a.y - pga.y)*nna.y); // normal distance from p1 to pg
              if(fabs(d2a)<EPS4 || fabs(d1a)<EPS4 || fabs(d3a)<EPS4 || fabs(nna.x)<EPS4 || fabs(nna.y)<EPS4){ // using \partial_T/\partial_n,
  
                  coord n= interface_normal3 (point, ff, hhh),p1,p2;
		  double alpha = plane_alpha (ff[], n);
		  double length = line_length_center (n, alpha, &p1);//2-dimension,p is the center point of the interface
		  double x1,y1;
		  x1 = x + p1.x*Delta, y1 = y + p1.y*Delta;  
		  double x2,y2;
		  plane_center(n,alpha,ff[],&p2);  // center of liquid
		  x2 = x + p2.x*Delta, y2 = y + p2.y*Delta;
		  coord dp1top2,nn;
		  dp1top2.x = x2 - x1;
		  dp1top2.y = y2 - y1;
		  nn.x = n.x/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
		  nn.y = n.y/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
		  double dd= dp1top2.x * nn.x + dp1top2.y * nn.y;
                              
                  Tlff[] = Tsat0 + phase1Tgrad[]*dd;
              }else{
                  Tlff[] = T[]*(d1a/d2a) - Tsat0*((d1a-d2a)/d2a);
              }
            }else if(ff[]<1.0){ 
                 Tlff[] = T[]; //COMMENT 20220812Tsat0;
                //Tlff[] = Tsat0;
            }else{ //pure liquid cell
                Tlff[] = T[];
           }
           
           if(ff[]<1.0-EPS){
              //for temperature of fration less than 0.5, using \partial_T/\partial_n, Tsat0;  becacue, when ff<0.5 we don't calculate Tgff in diffusion equation
              coord n= interface_normal3 (point, ff, hhh),p1,p2,p3;
	      double alpha = plane_alpha (ff[], n);
	      double length = line_length_center (n, alpha, &p1);//2-dimension,p is the center point of the interface
              double x1,y1;
	      x1 = x + p1.x*Delta, y1 = y + p1.y*Delta;  
              double x2,y2;
              plane_center(n,alpha,ff[],&p3);  // center of liquid
              p2.x = -ff[]*p3.x/(1.0-ff[]);
              p2.y = -ff[]*p3.y/(1.0-ff[]);
              x2 = x + p2.x*Delta, y2 = y + p2.y*Delta;
              coord dp1top2,nn;
              dp1top2.x = x2 - x1;
              dp1top2.y = y2 - y1;
              nn.x = n.x/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
              nn.y = n.y/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
              double dd= dp1top2.x * nn.x + dp1top2.y * nn.y;
                              
              Tgff[] = Tsat0 + phase0Tgrad[]*dd;
              //fprintf(fp57,"%d %d %g %g %g %g\n",point.i,point.j,x,y,Tlff[],Tgff[]);
             // fprintf(fp35,"g: i=%d j=%d ff=%g Tlff=%g Tgff=%g dd=%g nnx=%g nny=%g d12x=%g d12y=%g phase0Tgrad=%g \n",point.i,point.j, ff[], Tlff[], Tgff[], dd, nn.x, nn.y,dp1top2.x, dp1top2.y,phase0Tgrad[]);
            }else{
              Tgff[] = Tsat0; 
            }
      }
    
}else if(choice == 0){//choidce
    
	  for (scalar s in {Tlff,Tgff}){
		    if(s.inverse){
                        if(ff[]<0.5){
			   s[]=T[]*(1.0-ff[]);
                        }else{
                           s[]=Tsat0*ff[];
                        }
		    }else{
                        if(ff[]>=0.5){
			   s[]=T[]*ff[];
                        }else{
                           s[]=Tsat0*(1-ff[]);
                        }
		    }
         }

}
     
   }             // Linearly stratified
  //boundary ({T});

  //fclose(fp35);
 // fclose(fp56);
 // fclose(fp57);

if(choice != 0){
  foreach(){
	    for(scalar s in {Tlff,Tgff}){
		if(s.inverse){
		     s[]=Tgff[]*(1.0-ff[]);
		 }else{
	             s[]=Tlff[]*ff[];
		 }
	    }
   }
}

}

void  get_Tlff_Tgff_from_T_seperate(int phase_flag, scalar cc, scalar phaseTgrad){  // from T of cell center to bycenter of fluid.  This is from Tsat0 and \partial_T/partian_n
                       // another method is from T_sat T_liquid(for liquid cell) to approximate, for gas cell still using the Tsat0 and partial_T/partian_n
                       // scalar cc should be modified ff(according phase_flag)


double eps_new = 1e-10; //CCCCMIN;
double eps_new2 = 1e-10; //CCCCMIN;
  foreach(){
      
     if(cc[]>eps_new){  //liquid cell 
           if(cc[]<1.0-eps_new){  //from center temperature to fraction temperature , 
            
              ////for temperature of fration greater than 0.5, using T[], Tsat0
              coord na= interface_normal (point, cc),p1a;
	             double alphaa = plane_alpha (cc[], na);
             
              double lengtha = line_length_center (na, alphaa, &p1a);//2-dimension,p is centriod of the interface
              ////p1a.x = p1a.x*Delta, p1a.y= p1a.y*Delta;
              coord pla;
              plane_center(na,alphaa,cc[],&pla); //center of liquid
              ////pla.x = pla.x*Delta, pla.y= pla.y*Delta;
              coord pga;
              pga.x = -cc[]*pla.x/(1.0-cc[]);
              pga.y = -cc[]*pla.y/(1.0-cc[]);
              coord nna;
              nna.x = na.x/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              nna.y = na.y/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              double d1a,d2a,d3a;
              d2a = (p1a.x - 0.0)*nna.x +  (p1a.y - 0.0)*nna.y; // directional distance, toward oppsite normal of interface
              d1a = fabs((p1a.x - pla.x)*nna.x + (p1a.y - pla.y)*nna.y); // normal distance from p1 to pl
              d3a = fabs((p1a.x - pga.x)*nna.x + (p1a.y - pga.y)*nna.y); // normal distance from p1 to pg
             // if((fabs(d2a)<EPS4 || fabs(d1a)<EPS4 || fabs(d3a)<EPS4 || fabs(nna.x)<EPS4 || fabs(nna.y)<EPS4)){ // using \partial_T/\partial_n,
              if((fabs(d2a)<eps_new2)){// || fabs(d1a)<eps_new2) // || fabs(d3a)<eps_new2) // || fabs(nna.x)<eps_new2 || fabs(nna.y)<eps_new2)){
                    coord n= interface_normal3 (point, cc, hhh),p1,p2;
                    double alpha = plane_alpha (cc[], n);
                    double length = line_length_center (n, alpha, &p1);//2-dimension,p is the center point of the interface
                    double x1,y1;
                    x1 = x + p1.x*Delta, y1 = y + p1.y*Delta;  
                    double x2,y2;
                    plane_center(n,alpha,cc[],&p2);  // center of liquid
                    x2 = x + p2.x*Delta, y2 = y + p2.y*Delta;
                    coord dp1top2,nn;
                    dp1top2.x = x2 - x1;
                    dp1top2.y = y2 - y1;
                    nn.x = n.x/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
                    nn.y = n.y/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
                    double dd= dp1top2.x * nn.x + dp1top2.y * nn.y;
                    if(phase_flag == 1){          
                        Tlff[] = Tsat0 + phaseTgrad[]*dd;
                    }else{
                         Tgff[] = Tsat0 + phaseTgrad[]*dd;
                    }
              }else{
                    if(phase_flag == 1){
                        Tlff[] = T[]*(d1a/d2a) - Tsat0*((d1a-d2a)/d2a);
                    }else{
                         Tgff[] = T[]*(d1a/d2a) - Tsat0*((d1a-d2a)/d2a);
                    }
              }
            }else if(cc[]<1.0){ 
                // Tlff[] = T[]; //COMMENT 20220812Tsat0; ?????????????????????
                if(phase_flag == 1){
                     //Tlff[] = Tsat0;
                      Tlff[] = T[];
                }else{
                     //Tgff[] = Tsat0;
                      Tgff[] = T[];
                }
            }else{ //pure cell
                if(phase_flag == 1){
                    Tlff[] = T[];
                }else{
                    Tgff[] = T[];
                }
           }
      }else{ // full other liquid or small fraction liquid
            if(phase_flag == 1){
                  Tlff[] = Tsat0;
            }else{
                  Tgff[] = Tsat0;
            }

      }
    
     
   }             // Linearly stratified

   if(phase_flag == 1){
      foreach(){
          Tlff[] = T[] * cc[];
      }
   }else{
      foreach(){
          Tgff[] = T[] * cc[];
      }
   }

}



void  get_Tlff_Tgff_from_T2(){  // from T of cell center to bycenter of fluid.  This is from Tsat0 and \partial_T/partian_n
                       // another method is from T_sat T_liquid(for liquid cell) to approximate, for gas cell still using the Tsat0 and partial_T/partian_n

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// in diffusion step, Temperature in mixed cell is temperature in centerpoint temperature.
// in vof step, Temperature becomes temperature in bycenter(zhongxin) temperature. 
// This difference may be an improvement for the temperature diffusion equation.
// in ff>0.5 Tlff = temperature of the cell center point, while Tgff = Temprature of the bycenter(zhongxin)


//scalar phase0Tgrad[],phase1Tgrad[];  // 2022-0907-need to be deleted


    T_gradient_dirichlet(phase0Tgrad,phase1Tgrad);
//char name35[80];
//sprintf(name35,"Tlff_Tgff_mixed%g.dat", t);
//FILE * fp35 = fopen(name35,"a");

//char name56[80];
//sprintf(name56,"gas_cell_liquid_temp%g.dat", t);
//FILE * fp56 = fopen(name56,"a");
//fprintf(fp56,"#i j x y Tlff Tgff\n");
//char name57[80];
//sprintf(name57,"liquid_cell_gas_temp%g.dat", t);
//FILE * fp57 = fopen(name57,"a");
//fprintf(fp57,"#i j x y Tlff Tgff\n");


int choice=0;
  foreach(){
      
  if(choice == 1)   {
      if(ff[]< 0.5-EPS){ //gas cell
           if(ff[]>EPS){  //from center temperature to fraction temperature , 
             ////for temperature of fration greater than 0.5, using T[], Tsat0
              //coord na= interface_normal (point, ff),p1a;
	      //double alphaa = plane_alpha (ff[], na);
              ////p1a.x = p1a.x*Delta, p1a.y= p1a.y*Delta;
              //double lengtha = line_length_center (na, alphaa, &p1a);//2-dimension,p is centriod of the interface
              
              //coord pla;
              //plane_center(na,alphaa,ff[],&pla); //center of liquid
              ////pla.x = pla.x*Delta, p1a.y= pla.y*Delta;
              //coord pga;
              //pga.x = -ff[]*pla.x/(1.0-ff[]);
              //pga.y = -ff[]*pla.y/(1.0-ff[]);
              //coord nna;
              //nna.x = na.x/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              //nna.y = na.y/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              //double d1a,d2a,d3a;
              //d2a = (p1a.x - 0.0)*nna.x +  (p1a.y - 0.0)*nna.y; // directional distance, toward oppsite normal of interface
              //d1a = fabs((p1a.x - pla.x)*nna.x + (p1a.y - pla.y)*nna.y); // normal distance from p1 to pl
              //d3a = fabs((p1a.x - pga.x)*nna.x + (p1a.y - pga.y)*nna.y); // normal distance from p1 to pg
              //if(fabs(d2a)<EPS4  || fabs(d1a)<EPS4 || fabs(d3a)<EPS4 || fabs(nna.x)<EPS4 || fabs(nna.y)<EPS4){ // using \partial_T/\partial_n,   T[] from diffusion is not that accurate.
                  coord n= interface_normal3 (point, ff,hhh),p1,p2,p3;
	          double alpha = plane_alpha (ff[], n);
	          double length = line_length_center (n, alpha, &p1);//2-dimension,p is the center point of the interface
                  double x1,y1;
	          x1 = x + p1.x*Delta, y1 = y + p1.y*Delta;  
		  double x2,y2;
		  plane_center(n,alpha,ff[],&p3);  // center of liquid
		  p2.x = -ff[]*p3.x/(1.0-ff[]);
		  p2.y = -ff[]*p3.y/(1.0-ff[]);
		  x2 = x + p2.x*Delta, y2 = y + p2.y*Delta;
		  coord dp1top2,nn;
		  dp1top2.x = x2 - x1;
		  dp1top2.y = y2 - y1;
		  nn.x = n.x/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
		  nn.y = n.y/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
		  double dd= dp1top2.x * nn.x + dp1top2.y * nn.y;
                              
                  Tgff[] = Tsat0 + phase0Tgrad[]*dd;
              //}else{
              //    Tgff[] = Tsat0*((d2a + d3a)/(d2a)) - T[]*(d3a/(d2a));
              //}
            }else if(ff[]>0.0){ 
                Tgff[] = T[];//////// COMMENT 20220812 Tsat0;
            }else{//pure gas cell(ff<)
                Tgff[] = T[];
            }
 
           if(ff[]>EPS){  
              coord n= interface_normal3 (point, ff, hhh),p1,p2;
	      double alpha = plane_alpha (ff[], n);
	      double length = line_length_center (n, alpha, &p1);//2-dimension,p is the center point of the interface
              double x1,y1;
	      x1 = x + p1.x*Delta, y1 = y + p1.y*Delta;  
              double x2,y2;
              plane_center(n,alpha,ff[],&p2);  // center of liquid
              x2 = x + p2.x*Delta, y2 = y + p2.y*Delta;
              coord dp1top2,nn;
              dp1top2.x = x2 - x1;
              dp1top2.y = y2 - y1;
              nn.x = n.x/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
              nn.y = n.y/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
              double dd= dp1top2.x * nn.x + dp1top2.y * nn.y;
                              
              Tlff[] = Tsat0 + phase1Tgrad[]*dd;
              //fprintf(fp56,"%d %d %g %g %g %g\n",point.i,point.j,x,y,Tlff[],Tgff[]);
             // fprintf(fp35,"l: i=%d j=%d ff=%g Tlff=%g Tgff=%g dd=%g nnx=%g nny=%g d12x=%g d12y=%g phase1Tgrad=%g \n",point.i,point.j, ff[], Tlff[], Tgff[], dd, nn.x, nn.y,dp1top2.x, dp1top2.y,phase1Tgrad[]);
 
             }else{
                 Tlff[] = Tsat0;
             }  
             
      }else{  //liquid cell 
           if(ff[]<1.0-EPS){  //from center temperature to fraction temperature , 
            
              ////for temperature of fration greater than 0.5, using T[], Tsat0
              //coord na= interface_normal (point, ff),p1a;
	      //double alphaa = plane_alpha (ff[], na);
             
              //double lengtha = line_length_center (na, alphaa, &p1a);//2-dimension,p is centriod of the interface
              ////p1a.x = p1a.x*Delta, p1a.y= p1a.y*Delta;
              //coord pla;
              //plane_center(na,alphaa,ff[],&pla); //center of liquid
              ////pla.x = pla.x*Delta, pla.y= pla.y*Delta;
              //coord pga;
              //pga.x = -ff[]*pla.x/(1.0-ff[]);
              //pga.y = -ff[]*pla.y/(1.0-ff[]);
              //coord nna;
              //nna.x = na.x/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              //nna.y = na.y/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              //double d1a,d2a,d3a;
              //d2a = (p1a.x - 0.0)*nna.x +  (p1a.y - 0.0)*nna.y; // directional distance, toward oppsite normal of interface
              //d1a = fabs((p1a.x - pla.x)*nna.x + (p1a.y - pla.y)*nna.y); // normal distance from p1 to pl
              //d3a = fabs((p1a.x - pga.x)*nna.x + (p1a.y - pga.y)*nna.y); // normal distance from p1 to pg
              //if(fabs(d2a)<EPS4 || fabs(d1a)<EPS4 || fabs(d3a)<EPS4 || fabs(nna.x)<EPS4 || fabs(nna.y)<EPS4){ // using \partial_T/\partial_n,
  
                  coord n= interface_normal3 (point, ff, hhh),p1,p2;
		  double alpha = plane_alpha (ff[], n);
		  double length = line_length_center (n, alpha, &p1);//2-dimension,p is the center point of the interface
		  double x1,y1;
		  x1 = x + p1.x*Delta, y1 = y + p1.y*Delta;  
		  double x2,y2;
		  plane_center(n,alpha,ff[],&p2);  // center of liquid
		  x2 = x + p2.x*Delta, y2 = y + p2.y*Delta;
		  coord dp1top2,nn;
		  dp1top2.x = x2 - x1;
		  dp1top2.y = y2 - y1;
		  nn.x = n.x/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
		  nn.y = n.y/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
		  double dd= dp1top2.x * nn.x + dp1top2.y * nn.y;
                              
                  Tlff[] = Tsat0 + phase1Tgrad[]*dd;
              //}else{
              //    Tlff[] = T[]*(d1a/d2a) - Tsat0*((d1a-d2a)/d2a);
              //}
            }else if(ff[]<1.0){ 
                Tlff[] = T[]; //COMMENT 20220812Tsat0;
            }else{ //pure liquid cell
                Tlff[] = T[];
           }
           
           if(ff[]<1.0-EPS){
              //for temperature of fration less than 0.5, using \partial_T/\partial_n, Tsat0;  becacue, when ff<0.5 we don't calculate Tgff in diffusion equation
              coord n= interface_normal3 (point, ff, hhh),p1,p2,p3;
	      double alpha = plane_alpha (ff[], n);
	      double length = line_length_center (n, alpha, &p1);//2-dimension,p is the center point of the interface
              double x1,y1;
	      x1 = x + p1.x*Delta, y1 = y + p1.y*Delta;  
              double x2,y2;
              plane_center(n,alpha,ff[],&p3);  // center of liquid
              p2.x = -ff[]*p3.x/(1.0-ff[]);
              p2.y = -ff[]*p3.y/(1.0-ff[]);
              x2 = x + p2.x*Delta, y2 = y + p2.y*Delta;
              coord dp1top2,nn;
              dp1top2.x = x2 - x1;
              dp1top2.y = y2 - y1;
              nn.x = n.x/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
              nn.y = n.y/sqrt(pow(n.x,2.0)+pow(n.y,2.0));
              double dd= dp1top2.x * nn.x + dp1top2.y * nn.y;
                              
              Tgff[] = Tsat0 + phase0Tgrad[]*dd;
              //fprintf(fp57,"%d %d %g %g %g %g\n",point.i,point.j,x,y,Tlff[],Tgff[]);
             // fprintf(fp35,"g: i=%d j=%d ff=%g Tlff=%g Tgff=%g dd=%g nnx=%g nny=%g d12x=%g d12y=%g phase0Tgrad=%g \n",point.i,point.j, ff[], Tlff[], Tgff[], dd, nn.x, nn.y,dp1top2.x, dp1top2.y,phase0Tgrad[]);
            }else{
              Tgff[] = Tsat0; 
            }
      }
    
}else if(choice == 0){//choidce
    
	  for (scalar s in {Tlff,Tgff}){
		    if(s.inverse){
                        if(ff[]<0.5){
			   s[]=T[]*(1.0-ff[]);
                        }else{
                           s[]=Tsat0*ff[];
                        }
		    }else{
                        if(ff[]>=0.5){
			   s[]=T[]*ff[];
                        }else{
                           s[]=T[]*(1-ff[]);
                        }
		    }
         }

}
     
   }             // Linearly stratified
  //boundary ({T});

  //fclose(fp35);
 // fclose(fp56);
 // fclose(fp57);

if(choice != 0){
  foreach(){
	    for(scalar s in {Tlff,Tgff}){
		if(s.inverse){
		     s[]=Tgff[]*(1.0-ff[]);
		 }else{
	             s[]=Tlff[]*ff[];
		 }
	    }
   }
}

}

//#include "./myc2d2.h"
#include "./myc3d2.h"
/*
attribute {
  bool test3;
}
*/
void get_modphase01(){

    //if (ff.height.x.i){
     //     fprintf(stderr,"ff has height attribute!!!!!!\n");
     //}else{
    //      fprintf(stderr,"ff has no height attribute!!!!!!\n");
   //       heights(ff,ff.height);
     //      if (ff.height.x.i){
      //          fprintf(stderr,"then ff has\n");
      //     }
   //  }


    //first get n,aphsa,f and height; then restriction; for each face between phase=0 and phase=1, if there is no height
    // then calculate face-volume and normal,then get interface expretion and get the point
//scalar test1[];
    //ff5.test3 = true;

  //foreach_level(7){
    //if(child.x==1){
    //    fprintf(stderr,"childx exist !!!!!!!!\n");
    // }
  // }

//    face vector modphase1[];
//    face vector modphase0[];
    vector nn4[];
    scalar alpha4[];
    scalar ff4[];
  //  vector hhh[];
    foreach(){
       ff4[]=ff[];
    }   
    heights(ff,hhh);

//taken from restruction
foreach() {

    if (ff[] <= 0. || ff[] >= 1.) {
      alpha4[] = 0.;
      foreach_dimension()
	       nn4.x[] = 0.;
      }else{
      coord m = interface_normal3 (point, ff, hhh);
      foreach_dimension()
	        nn4.x[] = m.x;
      alpha4[] = plane_alpha (ff[], m);
    }
  }

#if TREE
  foreach_dimension()
    nn4.x.refine = nn4.x.prolongation = refine_injection;

  alpha4.n = nn4;
  alpha4.refine = alpha4.prolongation = alpha_refine;
#endif
  

  double lim_cut = 1e-12;
  restriction({ff4,nn4,alpha4});
  face vector complete1[];
  for(int l=0;l<=depth();l++){
       foreach_level(l){
          foreach_dimension(){
              complete1.x[] = 0;
              modphase1.x[] = HUGE;
              modphase0.x[] = HUGE;
              Point neib = neighborp(1);
              if(is_boundary(neighbor(1))){
                   foreach_neighbor(1){
                      if(point.i==neib.i && point.j==neib.j && point.k==neib.k){
                         modphase1.x[] = HUGE;
                         modphase0.x[] = HUGE;
                         complete1.x[] = 0;
                      }
                   }
              }
          }
       }
  }
  for(int l=0;l<=depth();l++){
     foreach_level(l){
       //if(complete1.x[]==0){  //commment on 20221210
         //coord complete={1,1};
         foreach_dimension(){
               if((ff4[-1]==0.0 && ff4[]==1.0) || (ff4[-1]==1.0 && ff4[]==0.0)){
                     modphase1.x[]=0.5;
                     modphase0.x[]=0.5;
                     //complete.x = 6;
                     complete1.x[] = 6;
                }else if( (is_phase(ff4[-1]) && !is_phase(ff4[])) || (!is_phase(ff4[-1]) && is_phase(ff4[]))){ // 
                   //complete.x=3;
                    complete1.x[]=3;
                   //bool complete=false;
		                 if((is_phase(ff4[-1]) && !is_phase(ff4[])) && orientation(hhh.x[])==0){ // get from height,(1 0)           
                        if(height(hhh.x[])>-1.0 && height(hhh.x[])<0.0){
		                      modphase0.x[]=fabs(height(hhh.x[]));
                          modphase1.x[]=1.0-modphase0.x[];
                         //complete.x=2;
                          complete1.x[]=2;
                         }
                       
		                 }else if((!is_phase(ff4[-1]) && is_phase(ff4[])) && orientation(hhh.x[])==1){ // 0 1
	                       if(height(hhh.x[])<0.0 && height(hhh.x[]>-1.0)){
                           modphase1.x[]=fabs(height(hhh.x[]));
                           modphase0.x[]=1.0-modphase1.x[];
                           //complete.x=2;
                           complete1.x[]=2;
                          }
		                 }
                  // if(complete.x==3){ // get from middle cell
                       
                   //}
                 if(complete1.x[]==3){ // get from middle cell
                     //bool complete=false;
		                 if((is_phase(ff4[-1]) && !is_phase(ff4[])) && orientation(hhh.x[-1])==0){ // get from height,(1 0)           
                        if(height(hhh.x[-1])>0.0 && height(hhh.x[-1])<1.0){
		                      modphase1.x[]=fabs(height(hhh.x[-1]));
                          modphase0.x[]=1.0-modphase1.x[];
                         //complete.x=2;
                          complete1.x[]=2;
                         }
                       
		                 }else if((!is_phase(ff4[-1]) && is_phase(ff4[])) && orientation(hhh.x[-1])==1){ // 0 1
	                       if(height(hhh.x[-1])>0.0 && height(hhh.x[-1]<1.0)){
                           modphase0.x[]=fabs(height(hhh.x[-1]));
                           modphase1.x[]=1.0-modphase0.x[];
                           //complete.x=2;
                           complete1.x[]=2;
                          }
		                 }
                   }
                 //      double lim_cut = 1e-13;
                  // if(complete.x==2 && (modphase1.x[]<lim_cut)){
                 if(complete1.x[]==2 && (modphase1.x[]<lim_cut)){
                             modphase1.x[] = lim_cut;
                             modphase0.x[] = 1.0 - modphase1.x[];
                             //complete.x = 5;
                             complete1.x[] = 5;
                 //  }else if(complete.x==2 && (modphase1.x[]>1.0-lim_cut)){
                    }else if(complete1.x[]==2 && (modphase1.x[]>1.0-lim_cut)){
                             modphase1.x[] = 1.0 - lim_cut;
                             modphase0.x[] = 1.0 - modphase1.x[];
                             //complete.x=5;
                             complete1.x[] = 5;
                   } 
                  
	           }else{
                  modphase1.x[] = 1.0;
                  modphase0.x[] = 1.0;
                  complete1.x[] = 7;
	          }
        }
        //  if(complete.x==3 && modphase1.x[]>HUGE/2.0){
          if(complete1.x[]==3 && modphase1.x[]>HUGE/2.0){
               double c_temp[3][3][3];
               for(int i0=-1;i0<2;i0++){
                  for(int j0=-1;j0<2;j0++){
                      for(int k0=-1;k0<2;k0++){
                            coord leftvolume,rightvolume;  //different from linear2.h, rifhtvolume is the right volume of -1; left volume if left volume of 0;

                            coord a_left,b_left;
                            coord a_right,b_right;
                            //coord n_temp= interface_normal (point, ff);
                            coord n_temp;
                            n_temp.x = nn4.x[i0,j0,k0], n_temp.y = nn4.y[i0,j0,k0], n_temp.z = nn4.z[i0,j0,k0];

                            coord n_temp0;
                            n_temp0.x = nn4.x[i0-1,j0,k0], n_temp0.y = nn4.y[i0-1,j0,k0], n_temp0.z = nn4.z[i0-1,j0,k0];
                            //double alpha_temp = plane_alpha (ff[], n_temp);
                            double alpha_temp = alpha4[i0,j0,k0];
                            double alpha_temp0 = alpha4[i0-1,j0,k0];
                            //foreach_dimension(){
                                
                                a_left=(coord){-0.5,-0.5,-0.5};
                                b_left=(coord){0.5,0.5,0.5};
                                
                                    
                                a_right=(coord){-0.5,-0.5,-0.5};
                                b_right=(coord){0.5,0.5,0.5};   
                              
                                b_left.x = 0.0;
                                leftvolume.x=rectangle_fraction(n_temp,alpha_temp,a_left,b_left)*(b_left.x-a_left.x);        
                                a_right.x = 0.0;
                                rightvolume.x=rectangle_fraction(n_temp0,alpha_temp0,a_right,b_right)*(b_right.x-a_right.x);  
                                                  
                                c_temp[i0+1][j0+1][k0+1] = rightvolume.x + leftvolume.x;
                          //}// foreach_dimension
                      } //for k0
                  } //for j0
                } //for i0

                //I am here
                coord n_middle = mycs3(c_temp);  //get normal from middle volume distribution
                double c_middle = c_temp[0+1][0+1][0+1];
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
          //if(complete.y==3 && modphase1.y[]>HUGE/2.0){
          if(complete1.y[]==3 && modphase1.y[]>HUGE/2.0){
               double c_temp[3][3][3];
               for(int i0=-1;i0<2;i0++){
                  for(int j0=-1;j0<2;j0++){
                      for(int k0=-1;k0<2;k0++){
                          coord leftvolume,rightvolume;  //different from linear2.h, rifhtvolume is the right volume of -1; left volume if left volume of 0;
                          coord a_left,b_left;
                          coord a_right,b_right;
                              //coord n_temp= interface_normal (point, ff);
                          coord n_temp;
                          n_temp.x = nn4.x[i0,j0,k0], n_temp.y = nn4.y[i0,j0,k0], n_temp.z = nn4.z[i0,j0,k0];
                          coord n_temp0;
                          n_temp0.x = nn4.x[i0,j0-1,k0], n_temp0.y = nn4.y[i0,j0-1,k0],n_temp0.z = nn4.z[i0,j0-1,k0];
                              //double alpha_temp = plane_alpha (ff[], n_temp);
                          double alpha_temp = alpha4[i0,j0,k0];
                          double alpha_temp0 = alpha4[i0,j0-1,k0];
                              //foreach_dimension(){
                              
                          a_left=(coord){-0.5,-0.5,-0.5};
                          b_left=(coord){0.5,0.5,0.5};
                          
                              
                          a_right=(coord){-0.5,-0.5,-0.5};
                          b_right=(coord){0.5,0.5,0.5};
                        
                          b_left.y = 0.0;
                          leftvolume.y=rectangle_fraction(n_temp,alpha_temp,a_left,b_left)*(b_left.y-a_left.y);        
                          a_right.y = 0.0;
                          rightvolume.y=rectangle_fraction(n_temp0,alpha_temp0,a_right,b_right)*(b_right.y-a_right.y);  
                                            
                          c_temp[i0+1][j0+1][k0+1] = rightvolume.y + leftvolume.y;
                            //}// foreach_dimension
                     } //for k0
                  } //for j0
                } //for i0
                coord n_middle = mycs3(c_temp);  //get normal from middle volume distribution
                double c_middle = c_temp[0+1][0+1][0+1];
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
         
        //  foreach_dimension(){
        //     if(is_boundary(neighbor(1))){
        //         modphase1.x[]=
        //     }
        //  }
          if(complete1.z[]==3 && modphase1.z[]>HUGE/2.0){
               double c_temp[3][3][3];
               for(int i0=-1;i0<2;i0++){
                  for(int j0=-1;j0<2;j0++){
                      for(int k0=-1;k0<2;k0++){
                          coord leftvolume,rightvolume;  //different from linear2.h, rifhtvolume is the right volume of -1; left volume if left volume of 0;
                          coord a_left,b_left;
                          coord a_right,b_right;
                              //coord n_temp= interface_normal (point, ff);
                          coord n_temp;
                          n_temp.x = nn4.x[i0,j0,k0], n_temp.y = nn4.y[i0,j0,k0], n_temp.z = nn4.z[i0,j0,k0];
                          coord n_temp0;
                          n_temp0.x = nn4.x[i0,j0,k0-1], n_temp0.y = nn4.y[i0,j0,k0-1],n_temp0.z = nn4.z[i0,j0,k0-1];
                              //double alpha_temp = plane_alpha (ff[], n_temp);
                          double alpha_temp = alpha4[i0,j0,k0];
                          double alpha_temp0 = alpha4[i0,j0,k0-1];
                              //foreach_dimension(){
                              
                          a_left=(coord){-0.5,-0.5,-0.5};
                          b_left=(coord){0.5,0.5,0.5};
                          
                              
                          a_right=(coord){-0.5,-0.5,-0.5};
                          b_right=(coord){0.5,0.5,0.5};
                        
                          b_left.z = 0.0;
                          leftvolume.z=rectangle_fraction(n_temp,alpha_temp,a_left,b_left)*(b_left.y-a_left.y);        
                          a_right.z = 0.0;
                          rightvolume.z=rectangle_fraction(n_temp0,alpha_temp0,a_right,b_right)*(b_right.y-a_right.y);  
                                            
                          c_temp[i0+1][j0+1][k0+1] = rightvolume.z + leftvolume.z;
                            //}// foreach_dimension
                     } //for k0
                  } //for j0
                } //for i0
                coord n_middle = mycs3(c_temp);  //get normal from middle volume distribution
                double c_middle = c_temp[0+1][0+1][0+1];
                double alpha_middle = plane_alpha(c_middle,n_middle);
                double test=0;;
                if(fabs(n_middle.z) > 1e-12){
	                   test =  alpha_middle/n_middle.z - (-0.5);
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

                modphase1.z[]=test;
                modphase0.z[]=1.0-modphase1.z[];
               // complete.y = 4; 
                complete1.z[] = 4;               
   //           
          }// complete.z==3
     
      ////} // //commment on 20221210
   }

/*
  // for halo cells   
  if(l>=1){
       foreach_halo(prolongation,l-1){
          foreach_dimension(){
              //Point neib = neighborp(-1);
              if(is_refined(neighbor(-1))){
                  foreach_child(){
                      if(child.x==-1){
                           // calculate c1
                           
                      }     
                  }

              }

          }
       }

  }
*/


if(l>=1){
 //     char name108[80];
//	    sprintf(name108,"halo_level%d_%g.dat",l, t);
//	    FILE * fp108 = fopen(name108,"w");
       foreach_halo(prolongation,l-1){
          double xp = x, yp = y, zp = z;
          //foreach_dimension(1){
              //Point neib = neighborp(-1);
          //    if(is_refined(neighbor(-1))){
                   
                  
                  coord childp;
                  foreach_child(){
                     foreach_dimension(){ 
                 //     if(child.x==-1){
                          // childp.x = x - 0.5*(Delta);
                          // childp.y = y;
                          // childp.x=xp-child.x*(0.5*Delta);
                          // childp.y=yp+child.y*(0.5*Delta);
 //                          fprintf(fp108,"%g %g %g %g\n",childp.x,childp.y,hhh.x[],ff[]); 
                          // fprintf(fp108,"%g %g %g\n",x,y,ff[]);       
                          modphase0.x[] = 1.0;
                          modphase1.x[] = 1.0;
                 //     } 
                     }    
                  }
                       
             //     xp = xp - 0.5*Delta;
             //     yp = y;
             //     fprintf(fp108,"%g %g %g %g %g %g\n",xp,yp,height(hhh.x[]),x,y,ff[]); 
        //      }

          //}
       }
   //   fclose(fp108);
  }



  
 
// 
     //to set right and top boundary condition for modphase0 and modphase1,but not sure foreach_face_level is right or not
     //set right and top boundary modphase1 and modphase0 to 1.0, so try to avoid to set bubble attached this boundary
    // for(int l=depth();l>=0;l--){  // foreach_face_level is not correct,because it loops several times
	     foreach_level(l){
                 foreach_dimension(){
                  Point neib2 = neighborp(1);
                 // if(modphase0.x[1]>HUGE/2.0 && is_active(cell) && is_boundary(neighbor(1))){
                    if(modphase0.x[1]>HUGE/2.0  && is_boundary(neighbor(1))){
                       foreach_neighbor(1){
                           // if(is_boundary(cell) && modphase0.x[]>HUGE/2.0 ){
                           if(is_boundary(cell) && point.i==neib2.i && point.j==neib2.j && point.k==neib2.k){
                               modphase0.x[] = 1.0;
                               modphase1.x[] = 1.0;
                           }
                       }//foreach_neibor
                   }
                 }//foreach_dimension
	     } //foreach_level
   //  } //for l

     // output modphase0 modphase1 and height to check
}
 
          //double lim_cut =1e-13;
     	  for(int l=0;l<=depth();l++){

  
	    //char name101[80];
	   // sprintf(name101,"height_level%d_%g.dat",l, t);
	    //FILE * fp101 = fopen(name101,"w");
	    foreach_level(l){
    /* 
		foreach_dimension(2){
		    int orien = orientation(hhh.x[]);
		    // if(hhh.x[] == nodata) {
		        
		    //}else{
             //            if((is_phase(ff4[]) && !is_phase(ff4[-1])) || (!is_phase(ff4[]) && is_phase(ff4[-1]))){
                           double x2=x-Delta/2.0, y2=y;
		        //   fprintf (fp101, "%g %g %g %d %g %d %g %g %g %g %g\n", x, y, ff4[], l,
// hhh.x[],orien,height(hhh.x[]),x2,y2,modphase0.x[],modphase1.x[]);
fprintf (fp101, "%g %g %g %d %g %g %g %g\n", x, y, ff4[], l,
x2,y2,modphase0.x[],modphase1.x[]);
             //      }
             //      }
		    //}
		 }
     */
                foreach_dimension(){
                     if(fabs(modphase0.x[])<lim_cut){
                          modphase0.x[]=lim_cut;
                          modphase1.x[]=1.0-modphase0.x[];
                     }else if(fabs(modphase1.x[])<lim_cut){
                          modphase1.x[]=lim_cut;
                          modphase0.x[]=1.0-modphase1.x[];
                     }
                }
	    }
	 //   fclose(fp101);
      // boundary_level({modphase0,modphase1},l);
	  }
     
  
  }



///for masstr[] source_pc[] 
//scalar masstr[];
extern scalar masstr;
extern scalar source_pc;
extern scalar vtr;
//extern double Tkg,Tkl,Trhog,Trhol,Tcpg,Tcpl,hfg;
extern double source_total,total_area;
//extern scalar source_ff1,source_ff2,ff_old1;
//#if TGRAD_LEON
#if 1
     #include "Tgrad-leon.h"
#endif 


void mass_transfer_rate(){
   //scalar masstr[];
   
   source_total= 0.0;
   total_area=0.0;
//#if TGRAD_LEON
#if 1
      Tgrad_leon(ff,phase0Tgrad,phase1Tgrad);
      printf("Trad_leon\n");
#else
       T_gradient_dirichlet(phase0Tgrad,phase1Tgrad);
       printf("Tgrad_classic");
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



   foreach(){
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
	//  foreach (reduction(+:area))
	    if (ff[] > EPS && ff[] < 1. - EPS) {
	      coord n = interface_normal (point, ff), p;
	      double alpha = plane_alpha (ff[], n);
               area = plane_area_center (n, alpha, &p);
	      //area += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
	    }

        //s_v(i,j,k) = mdot(i,j,k)*inv_drho*cell_area(i,j,k)/dh //area is a ratio, not a dimension number
        source_pc[] = masstr[]*(1.0/Trhog - 1.0/Trhol)*area/Delta;
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




//20221019-for moving interface by evaporate
//because interface cells have a same level, this will dramatically
//simplify the problem,since then we don't need to consider situation
//when interface is between resolution boundary
//important only allow interface with a same level
extern double dt;
void LevelSetShift2VOFChange(){
    
     scalar overshoot[];
     foreach(){
          deltac[] = 0.0;
     }

//char name37[80];
//sprintf(name37,"deltac1-%g.dat", t);
//FILE * fp37 = fopen(name37,"a");
//char name40[80];
//sprintf(name40,"dalpha-%g.dat", t);
//FILE * fp40 = fopen(name40,"a");
int flag_overshoot=0;
     foreach(reduction(max:flag_overshoot)){
              overshoot[]=0.0;
          // if(ff[]>EPS && ff[]<1.0-EPS){
	     if(ff[]>1e-6 && (ff[]<1-1e-6) && (css_test3_n[]>0.0)){
              coord n= interface_normal (point, ff);
	            double alpha_old = plane_alpha (ff[], n);
              double magn=sqrt(n.x*n.x + n.y*n.y + n.z*n.z); //2D-3D
              double magn1=fabs(n.x)+fabs(n.y)+fabs(n.z);
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
        position2.z[]=0;
        
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
              int xx=0,yy=0,zz=0;
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
               if(fabs(n.z)>normal_limit){
                  if(overshoot[]>0.0){  
                        zz = sign(n.z);        
                  }else{
                        zz = -sign(n.z);
                  }
               }
              if(xx==0 && yy==0 && zz==0){
                   fprintf (stderr, "WARNING: xx=0 && yy=0 && zz==0\n");
               }
 // for boiling case below should be take care!!!!!!!!!!!!!!!!!!!!!!!    
 //mpi point.i take care
      // if(point.i+xx>BGHOSTS2 && point.j+yy>BGHOSTS2 && point.i+xx<N+BGHOSTS2+1 && point.j+yy<N+BGHOSTS2+1){
      //   if(!is_boundary_box2(neighbor(xx,yy,zz))){
      if(!is_boundary_box2(neighbor(xx,yy,zz)) && (css_test3_n[])>0.0){
              if(fabs(intmask[xx,yy,zz]-3.0)>0.5){
		  // double magn=sqrt(n.x*n.x+n.y*n.y);
                   double alpha_old = plane_alpha (ff[xx,yy,zz], n); //using normal of interface and volume of accept cell
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
              deltac_temp[] =  ffnew - ff[xx,yy,zz];
              position2.x[]=xx;
              position2.y[]=yy;
              position2.z[]=zz;
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
          for(int k0=-1; k0<=1; k0=k0+1){
        //if neighbor has overshoot, judge if that overshoot belongs toyou
                  if(!(i0==0 && j0==0 && k0==0) ){
                      if(fabs(deltac_temp[i0,j0,k0])>0.0){
                                int i1=-round(position2.x[i0,j0,k0]);
                                int j1=-round(position2.y[i0,j0,k0]);
                                int k1=-round(position2.z[i0,j0,k0]);
                               // if(i0==i1  && j0==j1 && k0==k1) {		   
                                                                if(i0==i1  && j0==j1 && k0==k1 && (css_test3_n[i0,j0,k0]>0.0)) {	
                                      deltac[] = deltac[] + deltac_temp[i0,j0,k0];
                                              // flagi[i0,j0]=flagi[i0,j0]+1;
                                              flagi[] = flagi[] + 1;
                            //                  fprintf(fp38,"i=%d j=%d ff=%g i0=%d j0=%d i1=%d j1=%d pos2(i0,j0).x=%g pos2(i0,j0).y=%g deltac=%g deltac_temp(i0,j0)=%g \n",point.i,point.j,ff[],i0,j0,i1,j1,position2.x[i0,j0],position2.y[i0,j0],deltac[],deltac_temp[i0,j0]);
                                }
                      }
                 }
          } //k0
       } //j0
    } //i0
  }
}

 //fclose(fp38);    
//char name36[80];
//sprintf(name36,"deltac%g.dat", t);
//FILE * fp36 = fopen(name36,"a");

//char name39[80];
//sprintf(name39,"flagi%g.dat", t);
//FILE * fp39 = fopen(name39,"a");



    // // // // if(globali%outstep==0){
    // // // //     char name44[80];
    // // // //     sprintf(name44,"outfacets/error-dc%g.dat", t);
    // // // //     FILE * fp44 = fopen(name44,"a");
    // // // //     fprintf(fp44,"#i j deltac");

    // // // //     foreach(){

    // // // //       // if(fabs(intmask[]-3)>0.5){
    // // // //       //    int dd = round(7.89);
    // // // //         //    fprintf(fp36,"i=%d j=%d ff=%g dc_enter=%g dc=%g xx=%g yy=%g intmask=%g \n ",point.i,point.j,ff[],deltac_temp[],deltac[],position2.x[],position2.y[],intmask[]);
    // // // //       // }
    // // // //       //if(fabs(deltac_temp[])>0.0){
    // // // //       //    fprintf(fp39,"i=%d j=%d f=%g deltac_temp=%g flagi=%g \n",point.i,point.j,ff[],deltac_temp[],flagi[]);
    // // // //       //}

    // // // //       // if(fabs(deltac[])>0.0){
    // // // //           //fprintf(fp44,"%d %d %g\n",point.i,point.j,deltac[]);
    // // // //       // }

    // // // //       if(fabs(deltac[])>0.0){
    // // // //           if(level == level_interface && fabs(z-L0/2.0)<=Delta/2.0)
    // // // //               fprintf(fp44,"%g %g %g %g\n",x,y,ff[],deltac[]);
    // // // //       }
    // // // //     }
    // // // //     //fclose(fp36);
    // // // //     //fclose(fp39);    
    // // // //     fclose(fp44);
    // // // // }  

}

extern scalar topo_scalar;
extern scalar phase0Tg,phase1Tg;

void mass_record_out(char * name){
// if((globali%outstep)==0){
 // get_T_from_Tlff_Tgff3(1,ff); //20230110
	//char name59[80];
	//sprintf(name59,"mass_record-%g.dat", t);
	FILE * fp59 = fopen(name,"w");
       // fprintf(fp59,"#i j masstr source_pc\n");
        foreach(){
          
          //   if(ff[]<1.0-CCCCMIN && ff[]>CCCCMIN){
                coord n = interface_normal (point, ff), p;
	        double alpha = plane_alpha (ff[], n);
                double area = plane_area_center (n, alpha, &p);
                double cha =phase0Tgrad[] -phase1Tgrad[];

             // if((level == level_interface) && (fabs(z-L0/2.0)<=Delta/2.0)){
          //    if(fabs(topo_mask[])<=1){
            // if((fabs(y-L0/2.0)<=Delta/2.0 && y<L0/2.0) && (fabs(x-L0/2.0)<=Delta/2.0  && x<L0/2.0)){
              if(1==1){
         //     if((fabs(y-L0/2.0)<=Delta/2.0 && y<L0/2.0)){
              //if((level == level_interface) && (Tgff[]>3.0+0.1)){
            	   //fprintf(fp59,"%g %g %g %g %g %g %g %g %g %g %g %g %d %g %d\n",x,y,z,ff[],masstr[],source_pc[],Tgff[],Tlff[],area,phase0Tgrad[],phase1Tgrad[],cha,level,Delta,topo_mask[]);
                //  foreach_neighbor(1){
                //           fprintf(fp59,"       %g %g %g %g %g %g %g %g %g %g %g %d %g\n",x,y,ff[],masstr[],source_pc[],Tgff[],area,phase0Tgrad[],phase1Tgrad[],cha,level,Delta);
                //  }
                // fprintf(fp59,"%g %g %g %g %g %g %g\n",x,y,z,ff[],Tgff[],Tlff[],topo_mask[]);
               // if(fabs(topo_mask[])<=2){
              // //  //  fprintf(fp59,"%g %g %g %g %g %g %g %g %g %g %g %g %g\n",x,y,z,ff[],T[],topo_mask[],css_test[],css_test3[],css_test3_n[],phase0Tg[],phase1Tg[],DD.x[],fs_solid[]);
              fprintf(fp59,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",x,y,z,ff[],T[],topo_mask[],css_test[],css_test3[],css_test3_n[],phase0Tg[],phase1Tg[],fs_solid[],cm[],cm[]);
             
              } 
         //    }
             
             //if(fabs(source_pc[])>EPS)
             // fprintf(fp59,"%g %g %g %g %g %d %g\n",x,y,ff[],source_pc[],T[],level,Delta);
        }
        fclose(fp59);

   
//}

}

// void  mov_interface_dc(){

// double small_fraction_limit_T = 1e-12; //EPS;
//     // foreach(){
//     //     if(ff[]<small_fraction_limit_T && ff[]>0.0){ //this is just for evaporate, if cooling add >1e-6 situation
//     //        Tlff[]=0.0;
//     //        Tgff[]=Tgff[]/(1.0-ff[])*(1.0-0.0);
//     //        ff[]=0.0;
//     //     }else if(ff[]>1-small_fraction_limit_T && ff[]<1.0){
//     //        Tlff[]=Tlff[]/ff[]*1.0;
//     //        Tgff[]=0.0;
//     //        ff[]=1.0;
//     //     }
//     // }

// scalar f_old[];
// foreach(){
// //move dc + clamp(f[],0.0,1,0)
//         f_old[] = ff[];
//         ff[] = ff[] + deltac[];
	    
//    }
// foreach(){
// for(scalar s in {Tlff,Tgff}){
// 		if(s.inverse){
//                     if(1.0-ff[] > 0.0 && 1.0-f_old[] > 0.0){
//                         //s[] = (s[]/(1.0-f_old[])); //  s[]=Tgff[]*(1.0-ff[]);
//                         s[] = s[]; //  s[]=Tgff[]*(1.0-ff[]);
//                     }else if(1.0-ff[] <= 0.0 && 1.0-f_old[] > 0.0){
//                         s[] = Tsat00;
//                         ff[] = 1.0;
//                     }else if(1.0-ff[] > 0.0 && 1.0-f_old[] <= 0.0){
//                         s[] = Tsat00;
//                         ff[] = ff[];
//                     }else if(1.0-ff[] <= 0.0 && 1.0-f_old[] <= 0.0){
//                         s[] = Tsat00;
//                         ff[] = 1.0;
//                     }
                   
		   
// 		 }else{
//                      if(ff[] > 0.0 && f_old[] > 0.0){
//                         s[] = s[];
//                     }else if(ff[] <= 0.0 && f_old[] > 0.0){
//                         s[] = Tsat00;
//                         ff[] = 0.0;
//                     }else if(ff[] > 0.0 && f_old[] <= 0.0){
//                         s[] = Tsat00;
//                         ff[] = ff[];
//                     }else if(ff[] <= 0.0 && f_old[] <= 0.0){
// 			                 s[] = Tsat00;
//                         ff[] = 0.0;
//                     }
	            
// 		 }
// 	    }

// }

// // foreach(){
// // 	    for(scalar s in {Tlff,Tgff}){
// // 		if(s.inverse){
// // 		     s[]=s[]*(1.0-ff[]);
// // 		 }else{
// // 	             s[]=s[]*ff[];
// // 		 }
// // 	    }
// //    }

// }

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
                    if(ff[] > 0.0 && f_old[] <= 0.0 ){ //
                        val_tot = val_tot + Tsat00*ff[];
                        wei = wei + ff[];
                        val_tot = val_tot + T[]*(1.0-ff[]);
                        wei = wei + (1.0-ff[]);
                    }else if(1.0-ff[] > 0.0 && 1.0-f_old[] <= 0.0){
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

extern face vector ulf,uf;
extern scalar rho;
extern face vector alpha;
extern scalar ps;
extern face vector usf;
//extern scalar topo_mask;
extern double thick_init,L0;



void ulf_function(){

face vector facemask[];
    
    foreach_face(){
        facemask.x[]=1;
        //if(is_boundary_box2(point) || is_boundary_box2(neighbor(-1))){   
        //     facemask.x[]=0;
        //}
        if(is_boundary_box2(cell) || is_boundary_box2(neighbor(-1))){   
             facemask.x[]=0;
        }
     }
/*
if(periodic_top_bottom==1){
    foreach_face(y){
      if(is_boundary_box2(point) || is_boundary_box2(neighbor(0,-1))){   
             facemask.y[]=1;
        }
    }
}
*/
     foreach_face(){
	 ulf.x[] = 0.0;
     }
         
     
     //////////////////////////////////////////////////////////
     // ???assume that interface cell will not be coarsen or refine
     //////////////////////////////////////////////////////////
     foreach_face(){
         //vofnonmodule.f90 165-193
       if(facemask.x[]==1){
          if (level==level_interface){
              int correct1 =0;
              // if()
          //   if(topo_mask[]==0 || topo_mask[-1]==0 && (topo_mask[-1]<1 && topo_mask[]<1)){
                  //delete -2 situation = vof advection
              if((topo_mask[-1]<1 && topo_mask[]<1) && (topo_mask[]==0 || topo_mask[-1]==0)  ){
              // if(fabs(topo_mask[])<=2){ 
                      ulf.x[] = uf.x[] + usf.x[];
                }else{
                      ulf.x[] = uf.x[];
                }
          }else{
              ulf.x[] = uf.x[];
          }

        }
        if(fabs(fs.x[])<=0){
            ulf.x[] = 0.0; 
            usf.x[] = 0.0;
        }
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
/*
void delete_small_fraction(){
    double small_fraction_limit_T = 1e-12; //EPS;
    foreach(){
        if(ff[]<small_fraction_limit_T && ff[]>0.0){ //this is just for evaporate, if cooling add >1e-6 situation
           Tlff[]=0.0;
           Tgff[]=Tgff[]/(1.0-ff[])*(1.0-0.0);
           ff[]=0.0;
        }else if(ff[]>1-small_fraction_limit_T && ff[]<1.0){
           Tlff[]=Tlff[]/ff[]*1.0;
           Tgff[]=0.0;
           ff[]=1.0;
        }
    }
}
*/
