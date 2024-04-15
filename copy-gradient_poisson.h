//for two dimension
#include "fractions.h"

extern face vector fss_test;
extern scalar css_test;
extern vector hhh;
extern vector hhh_ff_oppo;
extern double total_area_diff,flux_total_diff;
extern bool flag_v0_nodata;
extern scalar ff;
extern double k0,k1;
extern double Tkg,Tkl;
extern double Tsat00;

#include "curvature.h" 



coord interface_normal5 (Point point, scalar c, vector h)
{
  coord n;
  if (!h.x.i || (n = height_normal (point, c, h)).x == nodata)
    n = mycs (point, c);
  return n;
}

coord facet_normal2 (Point point, scalar c, face vector s,vector h)
{
  if (s.x.i >= 0) { // compute normal from face fractions
    coord n;
    double nn = 0.;
    foreach_dimension() {
      n.x = s.x[] - s.x[1];
      nn += fabs(n.x);
    }
    if (nn > 0.)
      foreach_dimension()
	n.x /= nn;
    else
      foreach_dimension()
	n.x = 1./dimension;
    return n;
  }
  return interface_normal5 (point, c, h);
}

coord facet_normal3 (Point point, scalar c, face vector s,vector h)
{

  return interface_normal5 (point, c, h);
  // if (s.x.i >= 0) { // compute normal from face fractions
  //   coord n;
  //   double nn = 0.;
  //   foreach_dimension() {
  //     n.x = s.x[] - s.x[1];
  //     nn += fabs(n.x);
  //   }
  //   if (nn > 0.)
  //     foreach_dimension()
	// n.x /= nn;
  //   else
  //     foreach_dimension()
	// n.x = 1./dimension;
  //   return n;
  // }
 
}


#if dimension == 2
      #define face_condition2(fss_test, css_test)						\
        (fss_test.x[i,j] > 0.5 && fss_test.y[i,j + (j < 0)] && fss_test.y[i-1,j + (j < 0)] &&	\
        css_test[i,j] && css_test[i-1,j])


      //I am here
      foreach_dimension()
      static inline double embed_face_gradient2_x (Point point, scalar a, int i)
      {
        int j = sign(fss_test.x[i,1] - fss_test.x[i,-1]);
        assert (css_test[i] && css_test[i-1]);
        if (face_condition2 (fss_test, css_test))
          return ((1. + fss_test.x[i])*(a[i] - a[i-1]) +
            (1. - fss_test.x[i])*(a[i,j] - a[i-1,j]))/(2.*Delta);
        return (a[i] - a[i-1])/Delta;
      }

#else //dimension ==3
      // foreach_dimension()
      // static inline coord embed_face_barycentre_css_z (Point point, int i)
      // foreach_dimension()
      // coord embed_face_barycentre_css_z (Point point, int i)
      foreach_dimension()
      static inline coord embed_face_barycentre_css_z (Point point, int i)
      {
        // Young's normal calculation
        coord n1 = {0};
        double nn = 0.;
        scalar f = fss_test.z;
        foreach_dimension(2) {
          n1.x = (f[-1,-1,i] + 2.*f[-1,0,i] + f[-1,1,i] -
            f[+1,-1,i] - 2.*f[+1,0,i] - f[+1,1,i]);
          nn += fabs(n1.x);
        }
        if (!nn)
          return (coord){0.,0.,0.};
        foreach_dimension(2)
          n1.x /= nn;
        // Position `p` of the face barycentre
        coord n, p1, p;
        ((double *)&n)[0] = n1.x, ((double *)&n)[1] = n1.y;
        double alpha = line_alpha (f[0,0,i], n);
        line_center (n, alpha, f[0,0,i], &p1);
        p.x = ((double *)&p1)[0], p.y = ((double *)&p1)[1], p.z = 0.;
        return p;
      }

      #define face_condition3(fss_test, css_test)						\
        (fss_test.x[i,j,k] > 0.5 && (fss_test.x[i,j,0] > 0.5 || fss_test.x[i,0,k] > 0.5) &&	\
        fss_test.y[i,j + (j < 0),0] && fss_test.y[i-1,j + (j < 0),0] &&			\
        fss_test.y[i,j + (j < 0),k] && fss_test.y[i-1,j + (j < 0),k] &&			\
        fss_test.z[i,0,k + (k < 0)] && fss_test.z[i-1,0,k + (k < 0)] &&			\
        fss_test.z[i,j,k + (k < 0)] && fss_test.z[i-1,j,k + (k < 0)] &&			\
        css_test[i-1,j,0] && css_test[i-1,0,k] && css_test[i-1,j,k] &&				\
        css_test[i,j,0] && css_test[i,0,k] && css_test[i,j,k])


        foreach_dimension()
        static inline double embed_face_gradient2_x (Point point, scalar a, int i)
        {
        // foreach_dimension()
        // double embed_face_gradient2_x (Point point, scalar a, int i)
        // {
// when run with mpi, code exit with assert(); here I forced the code to coninue by set the volue to 0;
#if 0
          assert (css_test[i] && css_test[i-1]);
          coord p = embed_face_barycentre_css_x (point, i);
          // Bilinear interpolation of the gradient (see Fig. 1 of Schwartz et al., 2006)
          int j = sign(p.y), k = sign(p.z);
          if (face_condition3(fss_test, css_test)) {
            p.y = fabs(p.y), p.z = fabs(p.z);
            return (((a[i,0,0] - a[i-1,0,0])*(1. - p.y) +
              (a[i,j,0] - a[i-1,j,0])*p.y)*(1. - p.z) + 
              ((a[i,0,k] - a[i-1,0,k])*(1. - p.y) +
              (a[i,j,k] - a[i-1,j,k])*p.y)*p.z)/Delta;
          }
          return (a[i] - a[i-1])/Delta;
#else
          if(css_test[i] && css_test[i-1]){
                coord p = embed_face_barycentre_css_x (point, i);
                // Bilinear interpolation of the gradient (see Fig. 1 of Schwartz et al., 2006)
                int j = sign(p.y), k = sign(p.z);
                if (face_condition3(fss_test, css_test)) {
                  p.y = fabs(p.y), p.z = fabs(p.z);
                  return (((a[i,0,0] - a[i-1,0,0])*(1. - p.y) +
                    (a[i,j,0] - a[i-1,j,0])*p.y)*(1. - p.z) + 
                    ((a[i,0,k] - a[i-1,0,k])*(1. - p.y) +
                    (a[i,j,k] - a[i-1,j,k])*p.y)*p.z)/Delta;
                }
                return (a[i] - a[i-1])/Delta;
          }else{
               //return (a[i] - a[i-1])/Delta;
               return 0.0;
          }
#endif
        }

#endif


#define quadratic2(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

//extern scalar css_test; //scalar cs[];
//extern face vector fss_test; //face vector fs[];


// foreach_dimension()
// double dirichlet_gradient3_x (Point point, scalar s, scalar css_test,
// 					   coord n, coord p, double bc,
// 					   double * coef , face vector fss_test)
foreach_dimension()
static inline double dirichlet_gradient3_x (Point point, scalar s, scalar css_test,
					   coord n, coord p, double bc,
					   double * coef , face vector fss_test)
{
  foreach_dimension()
    n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fss_test.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fss_test.x[i + (i < 0),j] && fss_test.y[i,j] && fss_test.y[i,j+1] &&
	  css_test[i,j-1] && css_test[i,j] && css_test[i,j+1])
	v[l] = quadratic2 (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fss_test.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fss_test.y[i,j,k+m] || !fss_test.y[i,j+1,k+m] ||
	    !fss_test.z[i,j+m,k] || !fss_test.z[i,j+m,k+1] ||
	    !css_test[i,j+m,k-1] || !css_test[i,j+m,k] || !css_test[i,j+m,k+1])
	  defined = false;
      if (defined)
	// bi-quadratic2 interpolation
	v[l] =
	  quadratic2 (z,
		     quadratic2 (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic2 (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic2 (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
	break;
    }
  if (v[0] == nodata) {
 //   if (1==1) {

    /**
    This is a degenerate case, we use the boundary value and the
    cell-center value to define the gradient. */
	  flag_v0_nodata = true;
    d[0] = max(1e-3, fabs(p.x/n.x));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }

  /**
  For non-degenerate cases, the gradient is obtained using either
  second- or third-order estimates. */
  
  *coef = 0.;
  if (v[1] != nodata) // third-order gradient
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  return (bc - v[0])/(d[0]*Delta); // second-order gradient
}

double dirichlet_gradient3 (Point point, scalar s, scalar css_test,
			   coord n, coord p, double bc, double * coef, face vector fss_test)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradient3_x (point, s, css_test, n, p, bc, coef, fss_test);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradient3_x (point, s, css_test, n, p, bc, coef, fss_test);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradient3_y (point, s, css_test, n, p, bc, coef, fss_test);
  return dirichlet_gradient3_z (point, s, css_test, n, p, bc, coef, fss_test);
#endif // dimension == 3
  return nodata;
}



double embed_flux2 (Point point, scalar s, face vector mu, double * val, double bc)
{

  /**
  If the cell does not contain a fragment of embedded boundary, the
  flux is zero. */
  
  *val = 0.;
  if (css_test[] >= 1. || css_test[] <= 0.)
    return 0.;

  /**
  If the boundary condition is homogeneous Neumann, the flux is
  zero. */

  bool dirichlet ;
  dirichlet = true;
  double grad;
 // double grad = s.boundary[embed] (point, point, s, &dirichlet);
 // if (!grad && !dirichlet)
 //   return 0.;

  /**
  We compute the normal, area and barycenter of the fragment of embedded
  boundary contained within the cell. */

  //coord n = facet_normal (point, css_test, fss_test), p;
  //coord n = facet_normal2(point,css_test,fss_test,hhh);
 
  //coord n = facet_normal3(point,css_test,fss_test,hhh);
  coord n = interface_normal(point,ff);//height_normal (point, ff, hhh); // remenber using ff not css_test, because they are different in the solid
  coord p;
  double alpha = plane_alpha (ff[], n);
  double area = plane_area_center (n, alpha, &p);

  /**
  If the boundary condition is Dirichlet, we need to compute the
  normal gradient. */

  double coef = 0.;
  if (dirichlet) {
    normalize (&n);
    // grad = dirichlet_gradient3 (point, s, css_test, n, p, grad, &coef,fss_test);
    grad = dirichlet_gradient3 (point, s, css_test, n, p, bc, &coef,fss_test);
  }

  /**
  We retrieve the (average) value of $\mu$ without the metric. */
  
  double mua = 0., fa = 0.;
  if((!flag_v0_nodata) && fabs(bc-Tsat00)<1e-3){ 
     mua = Tkl;//k1;//Tkl;
     total_area_diff += area;
     flux_total_diff += area * grad;
     *val = - mua*grad*area/Delta;
     return 0.0;
   }
  foreach_dimension() {
    mua += mu.x[] + mu.x[1];
    fa  += fss_test.x[] + fss_test.x[1];
  }

   total_area_diff += area;
   flux_total_diff += area * grad;


  *val = - mua/(fa + SEPS)*grad*area/Delta;
  return - mua/(fa + SEPS)*coef*area/Delta;
}



// foreach_dimension()
// double dirichlet_gradient3_css_test2_x (Point point, scalar s, scalar css_test2,
// 					   coord n, coord p, double bc,
// 					   double * coef , face vector fss_test2)
foreach_dimension()
static inline double dirichlet_gradient3_css_test2_x (Point point, scalar s, scalar css_test2,
					   coord n, coord p, double bc,
					   double * coef , face vector fss_test2)
{
  foreach_dimension()
    n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fss_test2.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fss_test2.x[i + (i < 0),j] && fss_test2.y[i,j] && fss_test2.y[i,j+1] &&
	  css_test2[i,j-1] && css_test2[i,j] && css_test2[i,j+1])
	v[l] = quadratic2 (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fss_test2.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fss_test2.y[i,j,k+m] || !fss_test2.y[i,j+1,k+m] ||
	    !fss_test2.z[i,j+m,k] || !fss_test2.z[i,j+m,k+1] ||
	    !css_test2[i,j+m,k-1] || !css_test2[i,j+m,k] || !css_test2[i,j+m,k+1])
	  defined = false;
      if (defined)
	// bi-quadratic2 interpolation
	v[l] =
	  quadratic2 (z,
		     quadratic2 (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic2 (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic2 (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
	break;
    }
  if (v[0] == nodata) {
 //   if (1==1) {

    /**
    This is a degenerate case, we use the boundary value and the
    cell-center value to define the gradient. */
	  flag_v0_nodata = true;
    d[0] = max(1e-3, fabs(p.x/n.x));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }

  /**
  For non-degenerate cases, the gradient is obtained using either
  second- or third-order estimates. */
  
  *coef = 0.;
  if (v[1] != nodata) // third-order gradient
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  return (bc - v[0])/(d[0]*Delta); // second-order gradient
}

double dirichlet_gradient3_css_test2 (Point point, scalar s, scalar css_test2,
			   coord n, coord p, double bc, double * coef, face vector fss_test2)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradient3_css_test2_x (point, s, css_test2, n, p, bc, coef, fss_test2);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradient3_css_test2_x (point, s, css_test2, n, p, bc, coef, fss_test2);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradient3_css_test2_y (point, s, css_test2, n, p, bc, coef, fss_test2);
  return dirichlet_gradient3_css_test2_z (point, s, css_test2, n, p, bc, coef, fss_test2);
#endif // dimension == 3
  return nodata;
}


double embed_flux2_css_test2 (Point point, scalar s, face vector mu, double * val, double bc)
{

  /**
  If the cell does not contain a fragment of embedded boundary, the
  flux is zero. */
  
  *val = 0.;
  if (css_test2[] >= 1. || css_test2[] <= 0.)
    return 0.;

  /**
  If the boundary condition is homogeneous Neumann, the flux is
  zero. */

  bool dirichlet ;
  dirichlet = true;
  double grad;
 // double grad = s.boundary[embed] (point, point, s, &dirichlet);
 // if (!grad && !dirichlet)
 //   return 0.;

  /**
  We compute the normal, area and barycenter of the fragment of embedded
  boundary contained within the cell. */

  //coord n = facet_normal (point, css_test2, fss_test2), p;
  //coord n = facet_normal2(point,css_test2,fss_test2,hhh);
  //coord n = facet_normal3(point,css_test2,fss_test2,hhh); //use ff or ff_oppo
  coord n = interface_normal(point, ff_oppo); //height_normal(point, ff_oppo, hhh_ff_oppo);
  coord p;
  double alpha = plane_alpha (ff_oppo[], n);
  double area = plane_area_center (n, alpha, &p);

  /**
  If the boundary condition is Dirichlet, we need to compute the
  normal gradient. */

  double coef = 0.;
  if (dirichlet) {
    normalize (&n);
    // grad = dirichlet_gradient3 (point, s, css_test2, n, p, grad, &coef,fss_test2);
    grad = dirichlet_gradient3_css_test2 (point, s, css_test2, n, p, bc, &coef,fss_test2);
  }

  /**
  We retrieve the (average) value of $\mu$ without the metric. */
  
  double mua = 0., fa = 0.;
  if((!flag_v0_nodata) && fabs(bc-Tsat00)<1e-3){
     mua = Tkg;//k0;//Tkg;
     total_area_diff += area;
     flux_total_diff += area * grad;
     *val = - mua*grad*area/Delta;
     return 0.0;
   }

  foreach_dimension() {
    mua += mu.x[] + mu.x[1];
    fa  += fss_test2.x[] + fss_test2.x[1];
  }
  
 
  total_area_diff += area;
  flux_total_diff += area * grad;


  *val = - mua/(fa + SEPS)*grad*area/Delta;
  return - mua/(fa + SEPS)*coef*area/Delta;;
}

double embed_flux3 (Point point, scalar s, face vector mu, double * val, double bc)
{

  /**
  If the cell does not contain a fragment of embedded boundary, the
  flux is zero. */
  
  *val = 0.;
  if (css_test[] >= 1. || css_test[] <= 0.)
    return 0.;

  /**
  If the boundary condition is homogeneous Neumann, the flux is
  zero. */

  bool dirichlet ;
  dirichlet = true;
  double grad;
 // double grad = s.boundary[embed] (point, point, s, &dirichlet);
 // if (!grad && !dirichlet)
 //   return 0.;

  /**
  We compute the normal, area and barycenter of the fragment of embedded
  boundary contained within the cell. */

  //coord n = facet_normal (point, css_test, fss_test), p;
  //coord n = facet_normal2(point,css_test,fss_test,hhh);
  coord n = facet_normal3(point,css_test,fss_test,hhh);
  coord p;
  double alpha = plane_alpha (css_test[], n);
  double area = plane_area_center (n, alpha, &p);
   coord pla;
  plane_center(n,alpha,css_test[],&pla); //center of liquid
  /**
  If the boundary condition is Dirichlet, we need to compute the
  normal gradient. */

  double coef = 0.;
  if (dirichlet) {
    normalize (&n);
    // grad = dirichlet_gradient3 (point, s, css_test, n, p, grad, &coef,fss_test);
       coord delta_p;
       delta_p.x = p.x - pla.x;
       delta_p.y = p.y - pla.y;
       double distancebc = fabs(n.x*delta_p.x + n.y*delta_p.y);
       distancebc =  max(1e-3, distancebc);
    // if(distancebc > 1e-4){
    //    //grad = (bc - s[])/(distancebc*Delta); //bc/(d[0]*Delta);
       grad = bc/(distancebc*Delta);
       coef = - 1./(distancebc*Delta);
       
    //    //coef????????????????????????
    // }else{
    //    grad = dirichlet_gradient3 (point, s, css_test, n, p, bc, &coef,fss_test);

    // }
  }

  /**
  We retrieve the (average) value of $\mu$ without the metric. */
  
  double mua = 0., fa = 0.;
  foreach_dimension() {
    mua += mu.x[] + mu.x[1];
    fa  += fss_test.x[] + fss_test.x[1];
  }
  *val = - mua/(fa + SEPS)*grad*area/Delta;
  return - mua/(fa + SEPS)*coef*area/Delta;;
}


#undef face_gradient2_x
#define face_gradient2_x(a,i)					\
  (fss_test.x[i] < 1. && fss_test.x[i] > 0. ?			\
   embed_face_gradient2_x (point, a, i) :			\
   (a[i] - a[i-1])/Delta)

#undef face_gradient2_y
#define face_gradient2_y(a,i)					\
  (fss_test.y[0,i] < 1. && fss_test.y[0,i] > 0. ?		\
   embed_face_gradient2_y (point, a, i) :			\
   (a[0,i] - a[0,i-1])/Delta)

#undef face_gradient2_z
#define face_gradient2_z(a,i)					\
  (fss_test.z[0,0,i] < 1. && fss_test.z[0,0,i] > 0. ?		\
   embed_face_gradient2_z (point, a, i) :			\
   (a[0,0,i] - a[0,0,i-1])/Delta)


  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

