// foreach_dimension()
// static inline double dirichlet_gradient3_css_test2_x (Point point, scalar s, scalar css_test2,
// 					   coord n, coord p, double bc,
// 					   double * coef , face vector fss_test2)
extern scalar css_test,css_test2;
extern face vector fss_test, fss_test2;

double above_value = 0.05;//0.01;

#define quadratic2(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

foreach_dimension()
static inline double dirichlet_gradient3_css_test2_x (Point point, scalar s, scalar css_test2,
					   coord n, coord p, double bc,
					   double * method , face vector fss_test2)
{
//   coord n,p;
//   foreach_dimension(){
//     n.x = - n1->x;
//     p.x =   p1->x;
//   }
foreach_dimension(){
    n.x = -n.x;
}
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fss_test2.x[(n.x > 0.)])
      defined = false;
  if (defined){
    int lplus = 1;
    if(css_test[]>0.6){
        // lplus = 0;
        lplus = 1;
    }
    for (int l = 0; l <= 1; l++) {
      int i = (l + lplus)*sign(n.x);
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
  }
  if (v[0] == nodata) {
 //   if (1==1) {

    /**
    This is a degenerate case, we use the boundary value and the
    cell-center value to define the gradient. */
	//   flag_v0_nodata = true;
    // d[0] = max(1e-3, fabs(p.x/n.x));
    // *coef = - 1./(d[0]*Delta);
    // return bc/(d[0]*Delta);
    // flag_v0_nodata = true;
    d[0] = max(above_value/2.0, fabs(p.x/n.x));
    double value;
    value = (bc - s[])/(d[0]*Delta);
    // *coef = - 1./(d[0]*Delta);
    *method = 1;
    return value;
  }

  /**
  For non-degenerate cases, the gradient is obtained using either
  second- or third-order estimates. */
  
//   *coef = 0.;
  *method = 0;
  if (v[1] != nodata){ // third-order gradient
    *method = 3;
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  }
  *method = 2;
  return (bc - v[0])/(d[0]*Delta); // second-order gradient
}

foreach_dimension()
static inline double dirichlet_gradient3_css_test_x (Point point, scalar s, scalar css_test,
					   coord n, coord p, double bc,
					   double * method , face vector fss_test)
{
//   coord n,p;
//   foreach_dimension(){
//     n.x = - n1->x;
//     p.x =   p1->x;
//   }
foreach_dimension(){
    n.x = -n.x;
}
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fss_test.x[(n.x > 0.)])
      defined = false;
  if (defined){
    int lplus = 1;
    if(css_test[]>0.6){
        // lplus = 0;
        lplus = 1;
    }
    for (int l = 0; l <= 1; l++) {
      int i = (l + lplus)*sign(n.x);
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
  }
  if (v[0] == nodata) {
 //   if (1==1) {

    /**
    This is a degenerate case, we use the boundary value and the
    cell-center value to define the gradient. */
	//   flag_v0_nodata = true;
    // d[0] = max(1e-3, fabs(p.x/n.x));
    // *coef = - 1./(d[0]*Delta);
    // return bc/(d[0]*Delta);
    // flag_v0_nodata = true;
    d[0] = max(above_value/2.0, fabs(p.x/n.x));
    double value;
    value = (bc - s[])/(d[0]*Delta);
    // *coef = - 1./(d[0]*Delta);
    *method = 1;
    return value;
  }

  /**
  For non-degenerate cases, the gradient is obtained using either
  second- or third-order estimates. */
  
//   *coef = 0.;
  *method = 0;
  if (v[1] != nodata){ // third-order gradient
    *method = 3;
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  }
  *method = 2;
  return (bc - v[0])/(d[0]*Delta); // second-order gradient
}



// double dirichlet_gradient3_css_test2 (Point point, scalar s, scalar css_test2,
// 			   coord n, coord p, double bc, double * coef, face vector fss_test2)
// {
double dirichlet_gradient3_css_test2 (Point point, scalar s, scalar css_test2,
			   coord* n1, coord* p1, double bc, double * method, face vector fss_test2)
{
coord n,p;
  foreach_dimension(){
    n.x = n1->x;
    p.x = p1->x;
  }
#if dimension == 2
  // foreach_dimension()
    if (fabs(n.x) >= fabs(n.y)){
      return dirichlet_gradient3_css_test2_x (point, s, css_test2, n, p, bc, method, fss_test2);
    }else{
      return dirichlet_gradient3_css_test2_y (point, s, css_test2, n, p, bc, method, fss_test2);
    }
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradient3_css_test2_x (point, s, css_test2, n, p, bc, method, fss_test2);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradient3_css_test2_y (point, s, css_test2, n, p, bc, method, fss_test2);
  return dirichlet_gradient3_css_test2_z (point, s, css_test2, n, p, bc, method, fss_test2);
#endif // dimension == 3
  return nodata;
}

double dirichlet_gradient3_css_test (Point point, scalar s, scalar css_test,
			   coord* n1, coord* p1, double bc, double * method, face vector fss_test)
{
  coord n,p;
  foreach_dimension(){
    n.x = n1->x;
    p.x = p1->x;
  }
#if dimension == 2
  // foreach_dimension()
    if (fabs(n.x) >= fabs(n.y)){
      return dirichlet_gradient3_css_test_x (point, s, css_test, n, p, bc, method, fss_test);
    }else{
      return dirichlet_gradient3_css_test_y (point, s, css_test, n, p, bc, method, fss_test);
    }
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradient3_css_test_x (point, s, css_test, n, p, bc, method, fss_test);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradient3_css_test_y (point, s, css_test, n, p, bc, method, fss_test);
  return dirichlet_gradient3_css_test_z (point, s, css_test, n, p, bc, method, fss_test);
#endif // dimension == 3
  return nodata;
}

double embed_3_order_Tgrad_gas (Point point, scalar s, scalar f, double * method, double bc)
{

  /**
  If the cell does not contain a fragment of embedded boundary, the
  flux is zero. */
  
  *method = 0.;
  if (css_test2[] >= 1. || css_test2[] <= 0.)
    return 0.;

  /**
  If the boundary condition is homogeneous Neumann, the flux is
  zero. */

  bool dirichlet ;
  dirichlet = true;
  double grad;

  // coord n = interface_normal(point, f); //height_normal(point, ff_oppo, hhh_ff_oppo);
  coord n = mycs(point,f);
  coord p;
  double alpha = plane_alpha (f[], n);
  double area = plane_area_center (n, alpha, &p);

  *method = 0.;
  if (dirichlet) {
    // double method2;
    normalize (&n);
// grad = dirichlet_gradient3 (point, s, css_test2, n, p, grad, &coef,fss_test2);
grad = dirichlet_gradient3_css_test2 (point,s,css_test2,&n,&p,bc,method,fss_test2);
    // *method = method2;
  }
  return grad;
}

double embed_3_order_Tgrad_liquid (Point point, scalar s, scalar f, double * method, double bc)
{

  /**
  If the cell does not contain a fragment of embedded boundary, the
  flux is zero. */
  
  *method = 0.;
  if (css_test[] >= 1. || css_test[] <= 0.)
    return 0.;

  /**
  If the boundary condition is homogeneous Neumann, the flux is
  zero. */

  bool dirichlet ;
  dirichlet = true;
  double grad;
 
  // coord n = interface_normal(point,f);//height_normal (point, ff, hhh); // remenber using ff not css_test, because they are different in the solid
  coord n = mycs(point,f);
  coord p;
  double alpha = plane_alpha (f[], n);
  double area = plane_area_center (n, alpha, &p);

  /**
  If the boundary condition is Dirichlet, we need to compute the
  normal gradient. */

  *method = 0.;
  if (dirichlet) {
    // double method2;
    normalize (&n);
    // grad = dirichlet_gradient3 (point, s, css_test, n, p, grad, &coef,fss_test);
grad = dirichlet_gradient3_css_test (point, s, css_test, &n, &p, bc, method ,fss_test);
    // *method = method2;
  }
  return grad;
}